
# Python Imports
import argparse
from typing import Any, Dict, List, Optional, Type, Union, Generator
from contextlib import contextmanager

# Python Imports (for testing)
import pytest
from io import StringIO
from contextlib import redirect_stdout

def hierarchical_action_factory(root_attr: str) -> Type[argparse.Action]:
    """
    Factory creating a custom argparse.Action to store values in a nested dictionary.

    Parameters
    ----------
    root_attr : str
        The attribute name in the Namespace where the dictionary will be stored.

    Returns
    -------
    Type[argparse.Action]
        A NestedAction class tailored to the specific root attribute.
    """
    class NestedAction(argparse.Action):
        def __call__(self, parser: argparse.ArgumentParser, namespace: argparse.Namespace, 
                     values: Any, option_string: Optional[str] = None) -> None:
            if not hasattr(namespace, root_attr):
                setattr(namespace, root_attr, {})
            # Ensure it's a dict if it was pre-populated with something else
            if not isinstance(getattr(namespace, root_attr), dict):
                setattr(namespace, root_attr, {})

            store = getattr(namespace, root_attr)
            parts = self.dest.split('.')
            for part in parts[:-1]:
                store = store.setdefault(part, {})
            store[parts[-1]] = values
    return NestedAction

class BaseHierarchicalParser:
    """
    A wrapper for argparse to handle hierarchical configurations for tools and rules.

    Parameters
    ----------
    root_attr : str, optional
        The attribute name in the argparse Namespace to store the nested dict, 
        by default "tools_cfg".
    config_dict : dict, optional
        A dictionary (e.g., from JSON/YAML) used for default lookup and 
        context validation, by default None.
    **kwargs : Any
        Additional arguments passed to the underlying argparse.ArgumentParser.
    """
    def __init__(self, root_attr: str = "tools_cfg", config_dict: Optional[Dict[str, Any]] = None, **kwargs: Any) -> None:
        self.parser: argparse.ArgumentParser = kwargs.get('parser') or argparse.ArgumentParser(**kwargs)
        self.root_attr: str = root_attr
        self.config_dict: Optional[Dict[str, Any]] = config_dict
        self._ctx_path: List[str] = [] 
        self._ctx_container: Union[argparse.ArgumentParser, argparse._ArgumentGroup] = self.parser
        self._registered_paths: List[List[str]] = []


    def _validate_context(self, nodes: List[str]) -> None:
        """Hook for validation in subclasses/mixins. Default: allow everything."""
        pass

    def _validate_param(self, nodes: List[str], param_name: str, code_default: Any = argparse.SUPPRESS) -> None:
        """
        Hook for validation in subclasses/mixins. Default: allow everything.
        
        Parameters
        ----------
        nodes : list of str
            Current hierarchy.
        param_name : str
            Cleaned parameter name.
        code_default : Any, optional
            The default value provided in the add_param call.
        """
        pass

    @contextmanager
    def context(self, *nodes: str, group_name: Optional[str] = None) -> Generator['BaseHierarchicalParser', None, None]:
        """
        Context manager to group parameters under a specific hierarchy.

        Parameters
        ----------
        *nodes : str
            The hierarchy levels (e.g., "rule_name", "tool_name").
        group_name : str, optional
            If provided, creates a named argument group in the CLI help output.

        Yields
        ------
        BaseHierarchicalParser
            The parser instance with the context applied.

        Raises
        ------
        ValueError
            If a config_dict is provided but the requested nodes do not exist in it.
        """
        #if self.config_dict and self._get_nested_config(list(nodes)) is None:
        #    raise ValueError(f"Context {nodes} not found in provided config dictionary.")

        # use hook; meaningfull in combination with StrictConfigMixin
        self._validate_context(list(nodes))
        
        if list(nodes) not in self._registered_paths:
            self._registered_paths.append(list(nodes))

        old_path, old_container = self._ctx_path, self._ctx_container
        self._ctx_path = list(nodes)

        if group_name:
            self._ctx_container = self.parser.add_argument_group(group_name)
        else:
            self._ctx_container = self.parser

        try:
            yield self
        finally:
            self._ctx_path, self._ctx_container = old_path, old_container

    def add_param(self, *flags: str, **kwargs: Any) -> None:
        """
        Add a parameter to the current context. 
        
        Automatically generates a hierarchical flag (e.g., --rule-tool-param) 
        and links it to the nested dictionary.

        Parameters
        ----------
        *flags : str
            CLI flags (e.g., "-m", "--min-length").
        **kwargs : Any
            Standard argparse add_argument parameters (type, help, etc.).
        """

        long_flag_raw = max(flags, key=len).lstrip('-')
        clean_name = long_flag_raw.replace('-', '_')

        # use hook; meaningfull in combination with StrictConfigMixin
        self._validate_param(self._ctx_path, clean_name, kwargs.get('default', argparse.SUPPRESS))
        
        hier_flag = "--" + "-".join(self._ctx_path + [long_flag_raw]).replace("_", "-")
        all_flags = list(flags)
        if hier_flag not in all_flags:
            all_flags.append(hier_flag)

        # Default waterfall: CLI > JSON > Code
        if 'default' not in kwargs and self.config_dict:
            section = self._get_nested_config(self._ctx_path)
            if section and clean_name in section:
                kwargs['default'] = section[clean_name]

        # In case a default IS provided and it is stated in config_dict too,
        # it is accepted (however discouraged) that
        # coded-in-config, presumed DRY settings are overruled in code:
        #   parser.add_param("--min-len", type=int, default=10)
        # This would bypass any default value being stated in config_dict

        # Use SUPPRESS so we can detect if a value was provided/stored
        kwargs.setdefault('default', argparse.SUPPRESS)
        
        hierarchical_dest = ".".join(self._ctx_path + [clean_name])
        kwargs['action'] = hierarchical_action_factory(self.root_attr)
        kwargs['dest'] = hierarchical_dest
        self._ctx_container.add_argument(*all_flags, **kwargs)

    def parse_args(self, args: Optional[List[str]] = None) -> argparse.Namespace:
        """
        Parse command line arguments and merge with configuration defaults.

        Parameters
        ----------
        args : list of str, optional
            The arguments to parse, by default sys.argv[1:].

        Returns
        -------
        argparse.Namespace
            The namespace containing both vanilla arguments and the nested 
            hierarchical dictionary.
        """

        # 1. Standard parse (only CLI flags trigger NestedAction here)
        namespace = self.parser.parse_args(args)
        
        if not hasattr(namespace, self.root_attr):
            setattr(namespace, self.root_attr, {})
        target_dict = getattr(namespace, self.root_attr)

        # 2. Ensure all registered contexts exist in the dict
        for path in self._registered_paths:
            curr = target_dict
            for node in path:
                curr = curr.setdefault(node, {})

        # 3. Manually trigger actions for params NOT on CLI (applying defaults)
        for action in self.parser._actions:
            if "NestedAction" in str(type(action)):
                parts = action.dest.split('.')
                curr = target_dict
                for p in parts[:-1]: curr = curr[p]
                
                # If key not in dict, CLI didn't touch it. Apply default.
                if parts[-1] not in curr:
                    val = action.default if action.default is not argparse.SUPPRESS else None
                    curr[parts[-1]] = val

        # 4. Inject remaining JSON values for active branches
        if self.config_dict:
            self._inject_json_defaults(self.config_dict, target_dict)
            
        return namespace

    def _get_nested_config(self, path_list: List[str]) -> Optional[Dict[str, Any]]:
        """
        Traverse the config_dict to find a specific sub-section.

        Parameters
        ----------
        path_list : list of str
            The keys to follow in the configuration dictionary.

        Returns
        -------
        dict or None
            The sub-section if found, otherwise None.
        """
        if not self.config_dict: return None
        val = self.config_dict
        for node in path_list:
            if isinstance(val, dict) and node in val:
                val = val[node]
            else: return None
        return val if isinstance(val, dict) else None


    def _inject_json_defaults(self, source: Dict[str, Any], target: Dict[str, Any]) -> None:
        """
        Deep merge source configuration into the target dictionary.

        Parameters
        ----------
        source : dict
            The reference configuration (e.g., from a JSON file).
        target : dict
            The dictionary being populated by the parser.
        """
        for key, value in source.items():
            if key in target:
                if isinstance(value, dict) and isinstance(target[key], dict):
                    # If target is empty, we can just update it with the whole source dict
                    if not target[key]:
                        target[key].update(value)
                    else:
                        self._inject_json_defaults(value, target[key])
                elif target[key] is None:
                    target[key] = value

class StrictConfigMixin:
    """Implements the actual validation logic for strict parsers."""
    def _validate_context(self, nodes: List[str]) -> None:
        """
        Check if the requested context path exists in the configuration.

        Parameters
        ----------
        nodes : list of str
            The hierarchy levels to validate (e.g., ["rule_name", "tool_name"]).

        Raises
        ------
        ValueError
            If config_dict is None or if the path specified by `nodes` 
            does not exist in the configuration dictionary.
        """
        if self.config_dict is None:
            raise ValueError("Strict mode requires a config_dict.")
        if self._get_nested_config(nodes) is None:
            raise ValueError(f"Strict Mode: Context {nodes} not found in config.")

    def _validate_param(self, nodes: List[str], param_name: str, code_default: Any = argparse.SUPPRESS) -> None:
        """
        Check if a specific parameter exists in the configuration for the 
        current context.

        Parameters
        ----------
        nodes : list of str
            The current active hierarchy levels.
        param_name : str
            The cleaned name of the parameter to validate (e.g., "min_length").
        code_default : Any, optional
            The default value provided in the add_param call.

        Raises
        ------
        ValueError
            If the parameter is not found within the specified context 
            in the configuration dictionary.
        """
        section = self._get_nested_config(nodes)
        if section is None or param_name not in section:
            raise ValueError(f"Strict Mode: Parameter '{param_name}' missing in config for {nodes}.")

        # disallow overwriting a default value defined in the provided (exterior) config
        json_default = section.get(param_name)
        if code_default is not argparse.SUPPRESS and code_default != json_default:
            raise ValueError(
                f"Strict Mode: Default for '{param_name}' in code ({code_default}) "
                f"differs from JSON ({json_default}). JSON must be the source of truth."
            )


# from typing import Literal
# SeqType = Literal['illumina', 'nanopore', 'iontorrent', 'pacbio']
#
# def validate_literals(func):
#     """
#     Decorator that runtime-validates all arguments annotated with Literal.
#     """
#     @functools.wraps(func)
#     def wrapper(*args, **kwargs):
#         sig = inspect.signature(func)
#         bound_args = sig.bind(*args, **kwargs)
#         
#         for name, value in bound_args.arguments.items():
#             annotation = sig.parameters[name].annotation
#             # Check if it's a Literal type
#             if hasattr(annotation, "__origin__") and annotation.__origin__ is Literal:
#                 valid_choices = get_args(annotation)
#                 if value.lower() not in valid_choices:
#                     raise ValueError(
#                         f"Invalid value '{value}' for argument '{name}'. "
#                         f"Choices are: {valid_choices}"
#                     )
#                 # Normalize to lowercase automatically
#                 bound_args.arguments[name] = value.lower()
#         
#         return func(*bound_args.args, **bound_args.kwargs)
#     return wrapper
#
#class SeqTypeSnakemakeParser(BaseHierarchicalParser):
#    @contextmanager
#    @validate_literals
#    def tool_context(self, seq_type: SeqType, rule: str, tool: str, group_name: Optional[str] = None):
#        # De context manager handelt de 3 niveaus (nodes) automatisch af
#        with self.context(seq_type, rule, tool, group_name=group_name):
#            yield self
#

class SnakemakeParser(BaseHierarchicalParser):
    """Parser optimized for Snakemake workflows using [rule][tool] hierarchy."""
    @contextmanager
    def tool_context(self, rule: str, tool: str, group_name: Optional[str] = None):
        with self.context(rule, tool, group_name=group_name):
            yield self

class GenericToolParser(BaseHierarchicalParser):
    """Parser optimized for standalone tool configurations using [tool] hierarchy."""
    @contextmanager
    def tool_context(self, tool: str, group_name: Optional[str] = None):
        with self.context(tool, group_name=group_name):
            yield self

class ConfigDictatedSnakemakeParser(StrictConfigMixin, SnakemakeParser):
    """Identical to SnakemakeParser but with Strict validation via Mixin."""
    pass

class ConfigDictatedGenericToolParser(StrictConfigMixin, GenericToolParser):
    """Identical to GenericToolParser but with Strict validation via Mixin."""
    pass



# unit testing

@pytest.fixture
def bio_config():
    return {
        "fastq-QC": {
            "fastp": {"min_len": 50, "complexity_threshold": 30},
            "multiqc": {"title": "Raw_Data", "cl_config": "config.yaml"}
        },
        "mapping": {
            "bwa": {"min_seed_len": 19, "reseed_trigger": 1.5},
            "samtools": {"min_mq": 20, "mapq_filter": 30}
        },
        "mapping-QC": {
            "qualimap": {"homopolymer_size": 3, "genome_gc_dist": True},
            "picard": {"optical_dist": 100, "scoring_strategy": "SUM_OF_BASE_QUALITIES"}
        }
    }


@pytest.fixture
def snakemake_bio_parser(bio_config):
    """
    Provides a pre-configured SnakemakeParser with 3 rules/6 tools.
    Acts as a working example of the definition logic.
    """
    smk = SnakemakeParser(config_dict=bio_config)
    
    with smk.tool_context("fastq-QC", "fastp"):
        smk.add_param("--complexity-threshold", type=int)
    
    with smk.tool_context("mapping", "bwa"):
        smk.add_param("--reseed-trigger", type=float)
    
    with smk.tool_context("mapping-QC", "picard"):
        smk.add_param("--optical-dist", type=int)
        
    return smk

def test_full_hierarchical_flags(snakemake_bio_parser):
    """Validate that the auto-generated --rule-tool-param flags work."""
    args = snakemake_bio_parser.parse_args([
        "--mapping-bwa-reseed-trigger", "2.0",
        "--mapping-QC-picard-optical-dist", "250"
    ])
    
    assert args.tools_cfg["mapping"]["bwa"]["reseed_trigger"] == 2.0
    assert args.tools_cfg["mapping-QC"]["picard"]["optical_dist"] == 250
    # Ensure default is preserved for the untouched tool
    assert args.tools_cfg["fastq-QC"]["fastp"]["complexity_threshold"] == 30

def test_human_friendly_short_flags(snakemake_bio_parser):
    """Validate that the developer-defined short flags work and map correctly."""
    args = snakemake_bio_parser.parse_args([
        "--reseed-trigger", "2.5",
        "--optical-dist", "200"
    ])

    assert args.tools_cfg["mapping"]["bwa"]["reseed_trigger"] == 2.5
    assert args.tools_cfg["mapping-QC"]["picard"]["optical_dist"] == 200
    # Ensure default is preserved for the untouched tool
    assert args.tools_cfg["fastq-QC"]["fastp"]["complexity_threshold"] == 30

def test_mixed_flags_priority(snakemake_bio_parser):
    """Validate that if both are provided, the last one on the CLI wins (standard argparse)."""
    args = snakemake_bio_parser.parse_args([
        "--reseed-trigger", "3.0", 
        "--mapping-bwa-reseed-trigger", "4.0"
    ])
    assert args.tools_cfg["mapping"]["bwa"]["reseed_trigger"] == 4.0


def test_snakemake_config_injection(bio_config):
    """ Validate that defaults provided in the config are pulled into argument's default """
    smk = SnakemakeParser(config_dict=bio_config)
    with smk.tool_context("mapping", "bwa"):
        smk.add_param("--min-seed-len", type=int)
    
    args = smk.parse_args([])
    # Success: Default 19 is pulled and branch 'mapping' is created
    assert args.tools_cfg["mapping"]["bwa"]["min_seed_len"] == 19

def test_snakemake_cli_override(bio_config):
    """ Validate that one can provide an argument at the CLI and it overrules the stated default """
    smk = SnakemakeParser(config_dict=bio_config)
    with smk.tool_context("fastq-QC", "fastp"):
        smk.add_param("--min-len", type=int)
    
    # Success: Uses the hierarchical flag auto-generated by the engine
    args = smk.parse_args(["--fastq-QC-fastp-min-len", "75"])
    assert args.tools_cfg["fastq-QC"]["fastp"]["min_len"] == 75

def test_generic_tool_no_config():
    """ Validate that without a config dict the HierarchicalParser(s) still function """
    gtp = GenericToolParser(root_attr="custom_out")
    with gtp.tool_context("bwa"):
        gtp.add_param("--min-seed-len", type=int, default=10)
    
    args = gtp.parse_args([])
    assert args.custom_out["bwa"]["min_seed_len"] == 10

def test_generic_tool_with_config(bio_config):
    """ Validate that GenericToolParser basic functionality works """
    tool_cfg = bio_config["mapping"] 
    gtp = GenericToolParser(config_dict=tool_cfg)
    with gtp.tool_context("samtools"):
        gtp.add_param("--mapq-filter", type=int)
    
    args = gtp.parse_args(["--samtools-mapq-filter", "40"])
    assert args.tools_cfg["samtools"]["mapq_filter"] == 40

def test_injection_of_undefined_params(bio_config):
    """
    Test injection of JSON params that were NOT defined via add_param.
    This hits the 'elif target[key] is None' or full-dict injection logic.
    """
    smk = SnakemakeParser(config_dict=bio_config)
    
    # We define ONLY the rule context, but NO parameters via add_param
    # We manually trigger a branch in the dict to force the injector to work
    with smk.tool_context("fastq-QC", "fastp"):
        pass 
    
    # Because we opened the context, 'fastq-QC.fastp' exists in tools_cfg (as empty dict or None)
    # The injector should now fill it with EVERYTHING from the bio_config
    args = smk.parse_args([])
    
    assert "fastp" in args.tools_cfg["fastq-QC"]
    assert args.tools_cfg["fastq-QC"]["fastp"]["complexity_threshold"] == 30
    assert args.tools_cfg["fastq-QC"]["fastp"]["min_len"] == 50


def test_rare_param_not_in_json_is_none(bio_config):
    """
    Scenario: A 'rare' param is defined in code but undefined (deliberatly missing!) in JSON presets.
    User does not provide it on CLI. It should end up as None.
    """
    smk = SnakemakeParser(config_dict=bio_config)
    
    with smk.tool_context("mapping", "bwa"):
        # This param exists in code, but NOT in bio_config['mapping']['bwa']
        smk.add_param("--experimental-mode", type=str)
        # This one DOES exist in JSON
        smk.add_param("--min-seed-len", type=int)

    args = smk.parse_args([])
    
    # Check that the JSON default was still pulled
    assert args.tools_cfg["mapping"]["bwa"]["min_seed_len"] == 19
    # Check that the rare param is None and didn't crash the injector
    assert args.tools_cfg["mapping"]["bwa"]["experimental_mode"] is None

def test_coverage_trigger_elif_none(bio_config: Dict[str, Any]) -> None:
    """
    Forces the 'elif target[key] is None' branch using the bio_config fixture.
    Scenario: 
    1. A parser is initialized WITHOUT a config_dict first.
    2. add_param is called; since there is no config yet, it sets default to SUPPRESS.
    3. The config_dict is attached AFTER the parameters are defined.
    4. parse_args initializes the dict key to None (due to SUPPRESS).
    5. The injector sees the None and the value in bio_config, hitting the elif.
    """
    # 1. Initialize without the config_dict
    smk = SnakemakeParser(root_attr="tools_cfg", config_dict=None)
    
    # 2. Define a parameter that IS in the bio_config
    with smk.tool_context("mapping", "bwa"):
        # This will result in action.default = SUPPRESS because config_dict is None
        smk.add_param("--min-seed-len", type=int)

    # 3. NOW we attach the bio_config before parsing
    smk.config_dict = bio_config

    # 4. Parse
    args = smk.parse_args([])

    # Verify: The 'elif target[key] is None' must have triggered 
    # to move 19 from bio_config into the tools_cfg dict.
    assert args.tools_cfg["mapping"]["bwa"]["min_seed_len"] == 19

def test_context_creates_argument_group():
    """
    Verify that providing a group_name actually creates an 
    argparse._ArgumentGroup with the correct title.
    """
    p = SnakemakeParser()
    group_title = "BWA Mapping Options"
    
    with p.tool_context("mapping", "bwa", group_name=group_title):
        p.add_param("--threads", type=int)

    # Check internal action-groups for occurrence of the group_title
    # Index 0 and 1 are 'positional arguments' and 'options' in vanilla argparse
    group_titles = [g.title for g in p.parser._action_groups]
    
    assert group_title in group_titles
    
    # Check if the argument is indeed IN this particular action group
    target_group = next(g for g in p.parser._action_groups if g.title == group_title)
    group_dests = [a.dest for a in target_group._group_actions]
    assert "mapping.bwa.threads" in group_dests

def test_standard_parser_allows_new_hierarchical_entries(bio_config):
    """Confirm that the standard parser not raises ValueError for unknown contexts."""
    smk = SnakemakeParser(config_dict=bio_config)
    # This should now PASS without error
    with smk.tool_context("ghost-rule", "bwa"):
        smk.add_param("--new-param", type=int)

def test_vanilla_args_remain_top_level():
    """ Validate that any general parameter stays untouched in Namespace """
    smk = SnakemakeParser()
    smk.parser.add_argument("--verbose", action="store_true")
    with smk.tool_context("mapping", "bwa"):
        smk.add_param("--min-seed-len", type=int, default=19)
    
    args = smk.parse_args(["--verbose", "--mapping-bwa-min-seed-len", "25"])
    assert args.verbose is True
    assert args.tools_cfg["mapping"]["bwa"]["min_seed_len"] == 25

# tests various (other) edge cases

def test_nested_context_restoration():
    """
    Edge-case: Nested contexts.
    Ensures that when exiting a nested context, the path and container 
    return to the parent state, not to None.
    """
    p = SnakemakeParser()
    with p.tool_context("ruleA", "toolX") as p_outer:
        p_outer.add_param("--param1", type=int)
        
        with p.tool_context("ruleB", "toolY") as p_inner:
            p_inner.add_param("--param2", type=int)
            assert p._ctx_path == ["ruleB", "toolY"]
            
        # CRUCIAL: After exiting inner, path should be back to ruleA/toolX
        assert p._ctx_path == ["ruleA", "toolX"]
        p_outer.add_param("--param3", type=int)

    args = p.parse_args(["--ruleA-toolX-param1", "1", "--ruleB-toolY-param2", "2", "--ruleA-toolX-param3", "3"])
    assert args.tools_cfg["ruleA"]["toolX"]["param3"] == 3

def test_default_priority_logic(bio_config):
    """
    Edge-case: Priority between add_param(default=X) and config_dict.
    The rule should be: Explicit code default > JSON default.
    Mind this only holds true for non-StrictConfigMixin classes.
    Compare with test_strict_parsers_default_clash_comparison() 
    """
    # JSON has min_len: 50
    smk = SnakemakeParser(config_dict=bio_config)
    with smk.tool_context("fastq-QC", "fastp"):
        # We explicitly set a different default in code
        smk.add_param("--min-len", type=int, default=99)
    
    args = smk.parse_args([])
    # Code default should win over JSON
    assert args.tools_cfg["fastq-QC"]["fastp"]["min_len"] == 99

def test_nested_action_clobbers_non_dict_preexistence():
    """
    Forces coverage of the line: if not isinstance(getattr(namespace, root_attr), dict)
    Scenario: The namespace already has 'tools_cfg' but it's a string, not a dict.
    """
    p = SnakemakeParser(root_attr="tools_cfg")
    # We pre-populate the namespace with a string instead of a dict
    ns = argparse.Namespace(tools_cfg="I-am-a-string-garbage")
    
    with p.tool_context("rule", "tool"):
        p.add_param("--arg", type=int)
    
    # When we parse, the NestedAction must overwrite the string with a dict
    # We pass our 'dirty' namespace to the parser
    p.parser.parse_args(["--rule-tool-arg", "42"], namespace=ns)
    
    assert isinstance(ns.tools_cfg, dict)
    assert ns.tools_cfg["rule"]["tool"]["arg"] == 42

def test_add_param_on_provided_parser():
    """
    Test that a BaseHierarchicalParser can be initialized with an existing parser,
    and this existing and new added parameters are correctly added to the parser
    """
    base_parser = argparse.ArgumentParser(prog="test_tool")
    base_parser.add_argument(
        '--existing-flag',
        action="store_true",
        help="parameter defined in exterior ArgumentParser"
    )
    bio_config = {"some": "config"}
    smk = SnakemakeParser(parser=base_parser, config_dict=bio_config)
    # add an hierarchical parameter
    rule, tool ="fastqQC", "fastp"
    arg = "--min-length"
    with smk.tool_context(rule,tool):
        smk.add_param(arg, type=int, default=99)

    # add a "classical" parameter too
    base_parser.add_argument(
        '--new-classic-flag',
        action="store_true",
        help="classical parameter defined in exterior ArgumentParser"
    )

    # Check help output and search for a new and an earlier parameter
    f = StringIO()
    with redirect_stdout(f):
        try:
            smk.parse_args(['--help'])
        except SystemExit:
            pass
    help_output = f.getvalue()
    assert arg in help_output, "Parameter addition on SnakemakeParser(parser=base_parser) failed"
    assert "--"+("-".join([rule,tool,arg.lstrip('-')])) in help_output, "Parameter addition on SnakemakeParser(parser=base_parser) failed"
    assert "--existing-flag" in help_output, "Parameter from base_parser disappeared in SnakemakeParser(parser=base_parser)"
    assert "--new-classic-flag" in help_output, "Parameter from base_parser disappeared in SnakemakeParser(parser=base_parser)"

# tests for the ConfigDictatedMixin classes

def test_strict_snakemake_context_validation(bio_config):
    """Confirm ConfigDictatedSnakemakeParser blocks unknown rules/tools."""
    strict_smk = ConfigDictatedSnakemakeParser(config_dict=bio_config)
    with pytest.raises(ValueError, match="Strict Mode: Context.*not found"):
        with strict_smk.tool_context("ghost-rule", "bwa"):
            pass  # pragma: no cover

def test_strict_snakemake_param_validation(bio_config):
    """Confirm ConfigDictatedSnakemakeParser blocks unknown parameters."""
    strict_smk = ConfigDictatedSnakemakeParser(config_dict=bio_config)
    
    with strict_smk.tool_context("mapping", "bwa"):
        # This param is NOT in bio_config['mapping']['bwa']
        with pytest.raises(ValueError, match="Strict Mode: Parameter 'unknown_param' missing"):
            strict_smk.add_param("--unknown-param", type=int)

def test_strict_generic_tool_validation(bio_config):
    """Confirm ConfigDictatedGenericToolParser blocks unknown tools."""
    # Flatten config for generic tool
    tool_cfg = bio_config["mapping"]
    strict_gtp = ConfigDictatedGenericToolParser(config_dict=tool_cfg)
    
    with pytest.raises(ValueError, match="not found in config"):
        with strict_gtp.tool_context("ghost-tool"):
            pass # pragma: no cover

def test_strict_parsers_require_config():
    """Confirm that strict parsers raise error if initialized without config_dict."""
    with pytest.raises(ValueError, match="Strict mode requires a config_dict"):
        strict_smk = ConfigDictatedSnakemakeParser(config_dict=None)
        with strict_smk.tool_context("any", "thing"):
            pass # pragma: no cover

def test_strict_parsers_default_clash_comparison(bio_config):
    """
    Comparison test between Standard and Strict parsers regarding default overrides.
    """
    # 1. Standard Parser: Developer override is ALLOWED
    gtp = GenericToolParser(config_dict=bio_config["mapping"])
    with gtp.tool_context("bwa"):
        # JSON says 19, code says 99. Code wins in standard mode.
        gtp.add_param("--min-seed-len", type=int, default=99)
    
    args = gtp.parse_args([])
    assert args.tools_cfg["bwa"]["min_seed_len"] == 99

    # 2. Strict Parser: Developer override raises ValueError
    strict_gtp = ConfigDictatedGenericToolParser(config_dict=bio_config["mapping"])
    with strict_gtp.tool_context("bwa"):
        # This must raise ValueError because 99 != 19
        with pytest.raises(ValueError, match="differs from JSON"):
            strict_gtp.add_param("--min-seed-len", type=int, default=99)
