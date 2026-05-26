from __future__ import annotations

"""
Apollo mapping, based on Juno template
Authors: Roxanne Wolthuis, Boas van der Putten
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium
            Surveillance (IDS), Bacteriologie (BPD)     
Date: 07-04-2023   
"""

# Python Imports
from pathlib import Path
import yaml
import argparse
import os
import sys
import re
import pandas as pd
from tabulate import tabulate
from dataclasses import dataclass, field
from typing import Union, Optional, List, TYPE_CHECKING
# TODO chore: deprecate version.py file ?
from version import __package_name__, __version__

# Python Imports (RIVM packages)
import apollo_reference
from apollo_reference.referencedata import (
    ProvidedSchema,
    validate_reference_dataset,
    get_identify_species_mmidx_relpath,
)
from apollo_reference.dataframes import (
    read_reference_species_df,
    read_reference_assembly_df,
    read_reference_species_and_assembly_df,
    FASTQ_REFERENCE_TSV
)
from juno_library import Pipeline
from rivm_idsbioinfo_argparse_utils.actions import DynamicHelpTopicAction, DynamicHelpTopicShowMarkDownAction
from rivm_idsbioinfo_argparse_utils.parsers.hierarchicalconfigparser import SnakemakeParser
from rivm_idsbioinfo_argparse_utils.shortcuts import ArgumentContainer, get_flags_from_dest
from rivm_idsbioinfo_argparse_utils.argumentlabels import label_expert_arg, label_notimplemented_arg, label_deprecated_arg
from rivm_idsbioinfo_argparse_utils.formats.fasta import (
    validate_fasta_file,
    are_no_accession_overlap_in_fasta,
    are_accessions_present_in_fasta,
)

# configuration. read from various files from apollo-reference and apollo-mapping tool_parameters.yaml
APOLLO_REFERENCE_CONFIG_YAML = os.path.join(apollo_reference.__path__[0],"config","referencedata.yaml")
TOOL_PARAMS_YAML = os.path.join(os.path.dirname(__file__),"config","tool_parameters.yaml")
REFERENCE_DATA_DIR = yaml.safe_load(open(APOLLO_REFERENCE_CONFIG_YAML))['reference_data_dir']
KRAKEN_DB_DIR = yaml.safe_load(open(TOOL_PARAMS_YAML))['identify_impurity_using_kraken']['identify_impurity']['kraken_db']

CONFIG = {
    'kraken_db': Path(KRAKEN_DB_DIR),
    'apollo_reference_db_dir': Path(REFERENCE_DATA_DIR),
}

# _____________________________________________________________________________________________________ #
# TODO: move elsewhere
# _____________________________________________________________________________________________________ #


def as_argparse_type(validator):
    """ convert any validator (with a single argument) into an argparse type conversion function """
    def wrapper(value):
        try:
            return validator(value)
        except Exception as e:
            # Translate any Exception into an argparse error
            raise argparse.ArgumentTypeError(str(e))
    return wrapper

import functools
import inspect

def validate_arg(validator, arg_name=None):
    """ """
    def decorator(func):
        sig = inspect.signature(func)
        # define target_name: the given one OR defaulting to the first parameter od the function
        target_name = arg_name
        if target_name is None:
            # Take first (true instance method) parameter
            vanilla_class_parameters = ('self', 'cls')
            params = list(sig.parameters.keys())
            target_name = params[1] if params[0] in vanilla_class_parameters else params[0]

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()
            # Check if argument is there and validate!
            if target_name in bound_args.arguments:
                validator(bound_args.arguments[target_name])

            return func(*args, **kwargs)

        return wrapper

    return decorator

# _____________________________________________________________________________________________________ #
# TODO: move elsewhere
# _____________________________________________________________________________________________________ #

import six
import collections
def comma_separated_list(x: str ) -> list:
    """ split argument by comma's """
    if isinstance(x, six.string_types):
        # TODO: empty string returns ['']
        #       are empty strings allowed ?
        #       do we want empty string to get converted into None ?
        #       e.g. is "a,b,c,,,,,d,e" a valid argument" ?
        return x.split(',')
    elif isinstance(x, collections.abc.Iterable):
        # function used to validate an argument outside of argparse
        return list(x)
    else:
        raise ValueError("provided argument not string-type (neither iterable)")

# _____________________________________________________________________________________________________ #
# TODO: move elsewhere
# _____________________________________________________________________________________________________ #

def validate_kraken_db(path_str:Union[str|Path]):
    """Checks if the provided path represents a kraken2 directory """
    path = Path(path_str)

    # 1. If path doesn't exist or is no directory, return None
    if not path.is_dir():
        return None

    # 2. Check for required files for a kraked database directory
    required_files = [
        "database100mers.kmer_distrib",
        "database150mers.kmer_distrib",
        "database200mers.kmer_distrib",
        "database250mers.kmer_distrib",
        "database300mers.kmer_distrib",
        "database50mers.kmer_distrib",
        "database75mers.kmer_distrib",
        "hash.k2d",
        "inspect.txt",
        "ktaxonomy.tsv",
        "opts.k2d",
        "seqid2taxid.map",
        "taxo.k2d",
        ]
    if not all((path / f).exists() for f in required_files):
        msg = "provided argument lacks blueprint of a kraken2 database directory [%s]" % path_str
        raise ValueError(msg)

    return path

# ---------------------------------------------------------------------------------------- #
# helper functions
# ---------------------------------------------------------------------------------------- #

def slugify_str_choice_type(value: str, decapitalized:bool=True) -> str:
    """ helper function to slugify provided choice for species, accession, clade or cladegroup """
    value = value.strip()
    if re.search("^(GCA|GCF)_\d+", value):
        return value
    # slugify and de-capitalize
    value = value.replace(" ", "_").replace(".", "")
    if decapitalized:
        value = value[0].lower() + value[1:]
    return value

def slugify_str_choice_type_maintainedcase(value: str) -> str:
    """ helper function to slugify provided choice clade, maintaining the exact case usage """
    return slugify_str_choice_type(value,False)


class DynamicHelpSpeciesAction(DynamicHelpTopicAction):
    def generate_content(self) -> str:
        """ show default reference-dataset dataframe upon --help-species """
        df = read_reference_species_df()
        ProvidedSchema.validate(df, lazy=True)
        df.drop(columns=['ignore','literature','recent_name'], inplace=True)
        # pitch in MT origin from reference assembly information
        dfasm = read_reference_assembly_df()
        df = pd.merge(df, dfasm[['reference','taxid','MT_source']], on='reference', how='left')
        df['taxid'] = df['taxid'].astype('Int64')
        df.drop(columns=['mitochondrion'], inplace=True)
        df.rename(columns={'MT_source':'mitochondrion'}, inplace=True)
        return tabulate(df, headers='keys', showindex=False, intfmt="d", tablefmt='psql')

class DynamicHelpAssemblydataAction(DynamicHelpTopicAction):
    def generate_content(self) -> str:
        """ show default reference-assembly-dataset dataframe upon --help-assemblydata """
        # TODO: refactor minor non-DRY lines when compared to DynamicHelpSpeciesAction
        df = read_reference_species_and_assembly_df()
        df.drop(columns=['mitochondrion','MT_num','MT_assembly'], inplace=True)
        df.rename(columns={'MT_source':'MT'}, inplace=True)
        return tabulate(df, headers='keys', showindex=False, intfmt="d", tablefmt='psql')

class DynamicHelpCladesAction(DynamicHelpTopicAction):
    def generate_content(self) -> str:
        """ show manipulated default reference-dataset dataframe targeted for clade-analyses upon --help-clades """
        # TODO: refactor minor non-DRY lines when compared to DynamicHelpSpeciesAction
        df = read_reference_species_and_assembly_df()
        df = df[~df.cladegroup.isnull()]
        df.drop(columns=['mitochondrion','MT_num','MT_assembly'], inplace=True)
        df.rename(columns={'MT_source':'MT'}, inplace=True)
        df.sort_values(['cladegroup','is_primary','clade'],ascending=[True,False,True],inplace=True)
        return tabulate(df, headers='keys', showindex=False, intfmt="d", tablefmt='psql')

class DynamicHelpAvailableFastqAction(DynamicHelpTopicAction):
    def generate_content(self) -> str:
        """ show manipulated dataframe for references with publicly available fastq data (for testing) """
        # TODO: refactor minor non-DRY lines when compared to DynamicHelpSpeciesAction
        df = read_reference_species_and_assembly_df()
        df.drop(columns=['mitochondrion','MT_num','MT_assembly'], inplace=True)
        df.rename(columns={'MT_source':'MT'}, inplace=True)
        # load fastq dataframe
        dffq = pd.read_csv(FASTQ_REFERENCE_TSV,sep='\t')
        df = pd.merge(df, dffq, left_on='reference', right_on='accession', how='left')
        df.drop(columns=['accession'], inplace=True)
        df = df[~df.SRR .isnull()]
        return tabulate(df, headers='keys', showindex=False, intfmt="d", tablefmt='psql')

class DynamicHelpCustomrefsAction(DynamicHelpTopicShowMarkDownAction):
    MARKDOWN_FILE = "CHANGELOG.md"
    #def generate_content(self) -> str:
    #    """ show textual examples of how to apply custom reference data as input upon --help-customrefs """
    #    #text = read_markdown_file(path_to_xxxx)
    #    raise NotImplementedError("explanation to be developed")



def _add_argument_help_species(container: ArgumentContainer, flags: List[str]) -> None:
    """ add DynamicHelpSpeciesAction argument and link it to its argparse --flag """
    container.add_argument(
        *flags,
        action=DynamicHelpSpeciesAction,
        dest="help_species",
        help=f"show supported reference species and accession overview in tabular format",
    )

def _add_argument_help_assemblies(container: ArgumentContainer, flags: List[str]) -> None:
    """ add DynamicHelpAssemblydataAction argument and link it to its argparse --flag """
    container.add_argument(
        *flags,
        action=DynamicHelpAssemblydataAction,
        dest="help_assemblies",
        help=f"show supported reference species, with emphasis on the actual reference assemblies",
    )

def _add_argument_help_clades(container: ArgumentContainer, flags: List[str]) -> None:
    """ add DynamicHelpCladesAction argument and link it to its argparse --flag """
    container.add_argument(
        *flags,
        action=DynamicHelpCladesAction,
        dest="help_clades",
        help=f"show overview in tabular format of species with a defined cladegroup / clade subdivisions",
    )


def main() -> None:
    apollo_mapping = ApolloMapping()
    apollo_mapping.run()


@dataclass
class ApolloMapping(Pipeline):
    pipeline_name: str = __package_name__
    pipeline_version: str = __version__
    input_type: str = "fastq"
    #species_options = ["candida_auris", "aspergillus_fumigatus"]

    # options to validate provided -s / -a / -clg / -c <argument>
    species_options: list = field(default_factory=list)
    accession_options: list = field(default_factory=list)
    cladegroup_options: list = field(default_factory=list)
    clade_options: list = field(default_factory=list)

    def _set_argument_choices(self) -> None:
        """ assign choices for argparse arguments from default reference-tsv """
        df = read_reference_species_df()
        species = df.species.to_list()
        # mind type=str.LOWER in argument
        self.species_options = species + [ s.lower() for s in species ] + [ "_".join(s.lower().split()) for s in species ]
        self.species_options = list(sorted(set(self.species_options)))
        accessions = df.reference.unique().tolist()
        # temporarily patch for "eigen data". Will disappear once pandera validation is 100% strict
        try:
            accessions.pop(accessions.index("eigen data"))
        except ValueError:
            pass
        self.accession_options = accessions
        self.accession_options = list(sorted(set(self.accession_options)))

        self.cladegroup_options = df.cladegroup.dropna().unique().tolist()
        # TODO: clade should become unique
        # TODO: clade should be set when cladegroup is set
        clades = df.clade.dropna().unique().tolist()
        self.clade_options = clades + [ "_".join(c.split()) for c in clades ]

    def _assign_juno_library_pipeline_input_arguments_to_named_group(self):
        """ reshuffle generic Pipeline arguments into named groups [1/2] """
        # TODO: propagate this to juno-library. For now, monkey-patch based on the assumption
        #       of parameter ordering how "I" see in "now" on "my machine"
        # I know _assign_juno_library_pipeline_*_group is not DRY, no issue since will be solved at higher level
        param_index = next((i for i, a in enumerate(self.parser._actions) if "exclusion_file" == a.dest), None)
        group_title = 'Input Arguments'
        _ = self.parser.add_argument_group(group_title)
        group_input = self.parser._action_groups.pop()
        self.parser._action_groups.insert(0, group_input)

        actions_to_move = self.parser._actions[0:param_index+1]
        for action in actions_to_move:
            # add to this new group
            group_input._group_actions.append(action)
            for group in self.parser._action_groups:
                if action in group._group_actions and group.title != group_title:
                    group._group_actions.remove(action)

    def _assign_juno_library_pipeline_execution_arguments_to_named_group(self):
        """ reshuffle generic Pipeline arguments into named groups [2/2] """
        # TODO: propagate this to juno-library. For now, monkey-patch based on the assumption
        #       of parameter ordering how "I" see in "now" on "my machine"
        # I know _assign_juno_library_pipeline_*_group is not DRY, no issue since will be solved at higher level
        param_index_s = next((i for i, a in enumerate(self.parser._actions) if "exclusion_file" == a.dest), None)
        param_index_e = next((i for i, a in enumerate(self.parser._actions) if "snakemake_args" == a.dest), None)
        group_title = 'Execution Arguments: where and how do you want the pipeline executed?'
        group_execution = self.parser.add_argument_group(group_title)

        actions_to_move = self.parser._actions[param_index_s+1:param_index_e+1]
        for action in actions_to_move:
            # add to this new group
            group_execution._group_actions.append(action)
            for group in self.parser._action_groups:
                if action in group._group_actions and group.title != group_title:
                    group._group_actions.remove(action)
            # physically move these actions
            self.parser._actions.remove(action)
            self.parser._actions.append(action)


    def _move_help_to_custom_group(self) -> None:
        """ final monkey-patch: move vanilla --help to the argument_group with *--help* somewhere in its name """
        help_action = next((a for a in self.parser._actions if isinstance(a, argparse._HelpAction)), None)
        target_group = next((g for g in self.parser._action_groups if "--help" in g.title.lower()), None)

        if help_action and target_group:
            for group in self.parser._action_groups:
                if help_action in group._group_actions:
                    group._group_actions.remove(help_action)
            # place in group
            target_group._group_actions.append(help_action)
            # Optional: Physically move param to end of usage-list
            #parser._actions.remove(help_action)
            #parser._actions.append(help_action)


    def _add_args_to_parser(self) -> None:
        super()._add_args_to_parser()

        self.parser.description = (
            "Apollo mapping pipelines for reference mapping analysis of fungal genomes."
        )

        #bio_config = READ_RULE_TOOL_PARAM_FROM_JSON_OR_YAML()
        bio_config = {
            'map-clean-reads': {
                'picard': {
                    'ml': 50,
                }
            }
        }
        #bio_config = {}
        # prepare rule-tool-param specific overrides
        # after args = parser.parse_args(), access at args.workflowparams['map-clean-reads']['picard']['ml']
        self.smkp = SnakemakeParser(parser=self.parser,config_dict=bio_config,root_attr="workflowparams")
        self.smkp = SnakemakeParser(parser=self.parser,root_attr="workflowparams")
        #self.smkp = ConfigDictatedSnakemakeParser(parser=self.parser,config_dict=bio_config)

        self._set_argument_choices()
        self._assign_juno_library_pipeline_input_arguments_to_named_group()
        self._add_arguments_reference_species_scope()
        #self._add_arguments_tool_params()
        self._add_arguments_tool_params_hierarchicalconfigbased()
        self._add_arguments_tool_skipping()
        self._add_arguments_misc_help()
        self._add_arguments_deprecated_or_not_placed_yet()
        self._assign_juno_library_pipeline_execution_arguments_to_named_group()
        self._move_help_to_custom_group()

    def _add_arguments_reference_species_scope(self) -> None:
        """ """

        _group = self.parser.add_argument_group('Arguments changing reference species scope')
        mutexcl_group = _group.add_mutually_exclusive_group()

        mutexcl_group.add_argument(
            "-s",
            "--species",
            type=slugify_str_choice_type,
            metavar="STR",
            help=f"Reference Species to use; performs reference selection but dictates mapping to choosen species",
            dest="forced_species",
            choices=self.species_options,
        )

        mutexcl_group.add_argument(
            "-a",
            "--accession",
            type=slugify_str_choice_type,
            metavar="STR",
            help=f"Reference Accession to use; performs reference selection but dictates mapping to choosen accession",
            dest="forced_accession",
            choices=self.accession_options,
        )

        mutexcl_group.add_argument(
            "--clg",
            "--cladegroup",
            type=str,
            metavar="STR",
            help=f"Cladegroup to use; performs reference selection but dictates mapping to corresponding master reference",
            dest="forced_cladegroup",
            choices=self.cladegroup_options,
        )

        mutexcl_group.add_argument(
            "-c",
            "--clade",
            type=slugify_str_choice_type_maintainedcase,
            metavar="STR",
            help=f"Reference Clade to use; performs reference selection but dictates mapping to choosen Clade (and corresponding master reference)",
            dest="forced_clade",
            choices=self.clade_options,
        )

        _add_argument_help_species(mutexcl_group,["--species-help"])
        _add_argument_help_assemblies(mutexcl_group,["--assemblies-help"])
        _add_argument_help_clades(mutexcl_group,["--clades-help"])

        mutexcl_group.add_argument(
            "--custom-reference",
            type=as_argparse_type(validate_fasta_file),
            metavar="FASTA",
            default=None,
            dest="custom_reference_fasta",
            help="Custom Reference genome to use; overrules the reference selection & directly jumps to this reference",
        )

        label_notimplemented_arg(mutexcl_group.add_argument)(
            "--extra-reference",
            type=as_argparse_type(validate_fasta_file),
            metavar="FASTA",
            dest="extra_reference_fasta",
            help="Extra Reference genome to use; adds the provided reference (just once) to the reference selection step",
        )

        label_notimplemented_arg(mutexcl_group.add_argument)(
            "--custom-reference-tsv",
            type=Path, # is_valid_tsv_file
            metavar="FILE",
            default=None,
            dest="custom_reference_tsv",
            help="Custom Reference genomes to use (in tsv format); overrules the default supported reference set",
        )

        # TODO: this one should rewite SPECIES_REFERENCE_TSV and ASSEMBLY_REFERENCE_TSV
        label_expert_arg(mutexcl_group.add_argument)(
            "--custom-reference-dataset",
            type=Path, # is_valid_reference_dataset
            metavar="[PATH]",
            default=None,
            dest="custom_reference_dataset",
            help="Custom Reference Dataset directory to use; overrules the default supported reference set."
        )

        _group = self.parser.add_argument_group('Extra optional arguments in combination with --custom-reference')
        extra_customref_mito_group = _group.add_mutually_exclusive_group(required=False)
        extra_customclade_mito_group = _group.add_mutually_exclusive_group(required=False)


        extra_customref_mito_group.add_argument(
            "--custom-reference-mitochondrion-fasta",
            type=as_argparse_type(validate_fasta_file),
            metavar="FASTA",
            default=None,
            dest="custom_reference_mitochondrion_fasta",
            help="Custom Mitochondrion genome connected to your custom reference genome",
        )

        extra_customref_mito_group.add_argument(
            "--custom-reference-mitochondrion-accessions",
            type=as_argparse_type(comma_separated_list),
            metavar="[acc,...]",
            default=None,
            dest="custom_reference_mitochondrion_accessions",
            help="Accession(s) occurring in your --custom-reference that represent the mitochondrion"
        )


        label_notimplemented_arg(_group.add_argument)(
            "--custom-clade",
            type=as_argparse_type(validate_fasta_file),
            metavar="FASTA",
            default=None,
            dest="custom_clade_fasta",
            help="Custom Clade genome to use alongside a custom reference; overrules the reference/clade selection",
        )

        label_notimplemented_arg(extra_customclade_mito_group.add_argument)(
            "--custom-clade-mitochondrion-fasta",
            type=as_argparse_type(validate_fasta_file),
            metavar="FASTA",
            default=None,
            dest="custom_clade_mitochondrion_fasta",
            help="Custom Mitochondrion genome connected to your custom clade genome",
        )

        label_notimplemented_arg(extra_customclade_mito_group.add_argument)(
            "--custom-clade-mitochondrion-accessions",
            type=as_argparse_type(comma_separated_list),
            metavar="[acc,...]",
            default=None,
            dest="custom_clade_mitochondrion_accessions",
            help="Accession(s) occurring in your --custom-clade that represent the mitochondrion"
        )



    def _add_arguments_tool_params(self) -> None:
        tool_argument_group = self.parser.add_argument_group('Arguments that overrule arguments for tools in rules')

        tool_argument_group.add_argument(
            "-mpt",
            "--mean-quality-threshold",
            type=int,
            metavar="INT",
            default=28,
            # note: waarom
            choices=range(20,41),
            help="Phred score to be used as threshold for cleaning (filtering) fastq files.",
        )
        tool_argument_group.add_argument(
            "-ws",
            "--window-size",
            type=int,
            metavar="INT",
            default=5,
            help="Window size to use for cleaning (filtering) fastq files.",
        )
        tool_argument_group.add_argument(
            "-ml",
            "--minimum-length",
            type=int,
            metavar="INT",
            default=50,
            help="Minimum length for fastq reads to be kept after trimming.",
        )

    def _add_arguments_tool_params_hierarchicalconfigbased(self) -> None:
        """ """
        with self.smkp.tool_context("ithinkmultiqc", "icannotfindthemptflag"):
            self.smkp.add_param(
                "-mpt",
                "--mean-quality-threshold",
                type=int,
                metavar="INT",
                default=28,
                help="Phred score to be used as threshold for cleaning (filtering) fastq files.",
            )

        with self.smkp.tool_context("map_clean_reads", "picard"):
            self.smkp.add_param(
                "-ws",
                "--window-size",
                type=int,
                metavar="INT",
                default=5,
                help="Window size to use for cleaning (filtering) fastq files.",
            )
            self.smkp.add_param(
                "-ml",
                "--minimum-length",
                type=int,
                metavar="INT",
                default=50,
                help="Minimum length for fastq reads to be kept after trimming.",
            )


    def _add_arguments_tool_skipping(self) -> None:

        skip_tool_group = self.parser.add_argument_group('Arguments that disable particular tools in rules')

        #_group = skip_tool_group.add_mutually_exclusive_group(required=False)
        skip_tool_group.add_argument(
            "--db-dir",
            type=as_argparse_type(validate_kraken_db),
            default=CONFIG['kraken_db'],
            metavar="DIR",
            help="Kraken2 database directory (should include fungi!); use --skip-kraken to omit it",
        )

        skip_tool_group.add_argument(
            "--skip-kraken",
            action="store_true",
            default=False,
            help="omit the Kraken analyses (which quantifies possible sample impurity)"
        )

        label_notimplemented_arg(skip_tool_group.add_argument)(
            "--skip-reference-selection",
            action="store_true",
            default=False,
            help="skip the reference selection analyses (requires one of -s/-a/-c/--clg)"
        )

        label_notimplemented_arg(skip_tool_group.add_argument)(
            "--skip-multiclade-strain-mapping",
            action="store_true",
            default=False,
            help="skip, if applicable, the strain-specific mapping+variant calling on the designated strain (from a multiclade group)"
        )


        label_notimplemented_arg(skip_tool_group.add_argument)(
            "--trigger-assembly",
            action="store_true",
            default=False,
            help="trigger the Spades de-novo assembly of your sample(s) in the apollo-assembly pipeline"
        )

        argument_group = self.parser.add_argument_group('Particular use-case options')
        ISO_group = argument_group.add_mutually_exclusive_group(required=False)

        ISO_group.add_argument(
            "--ISO",
            action="store_true",
            dest="ISO",
            default=False,
            help="Fallback to ISO-certified part(s) of the pipeline"
        )

        ISO_group.add_argument(
            "--ISO-Cauris",
            action="store_true",
            dest="ISO_Cauris",
            default=False,
            help="Fallback to ISO-certified part(s) of the pipeline for C.auris (and disallow non-compatible arguments)"
        )


    def _add_arguments_misc_help(self) -> None:

        misc_help_group = self.parser.add_argument_group('Miscellaneous --help-<topic> that exhibits explanation on extra arguments')

        misc_help_group.add_argument(
            "--help-toolargs",
            action='version',
            dest="help_toolargs",
            version="xxxxxxxxxxxxxxxxxxxxxxxxxx",
            help=f"show overview in the various tool's arguments that can be adjusted in this pipeline",
        )
        misc_help_group.add_argument(
            "--help-customrefs",
            action=DynamicHelpCustomrefsAction,
            dest="help_customrefs",
            help="show some brief examples on how to parameterize input for custom reference data",
        )
        # realize (e.g.) --help-species and --species-help exist
        _add_argument_help_species(misc_help_group,["--help-species"])
        _add_argument_help_assemblies(misc_help_group,["--help-assemblies"])
        _add_argument_help_clades(misc_help_group,["--help-clades"])
        misc_help_group.add_argument(
            "--help-fastqdata",
            action=DynamicHelpAvailableFastqAction,
            dest="help_fastqdata",
            help="show reference data for which (public) fastq data is available",
        )



    def _add_arguments_deprecated_or_not_placed_yet(self) -> None:
        """ """
        etc_group = self.parser.add_argument_group('Miscellaneous: deprecated, unplaced ... the final loose ends')
        label_deprecated_arg(etc_group.add_argument)(
            "--reference",
            type=Path,
            metavar="FILE",
            dest="custom_reference",
            help="Reference genome to use default is chosen based on species argument, defaults per species can be found in: /mnt/db/apollo/mapping/[species]",
            required=False,
        )


    def _get_mmidx_path(self) -> str:
        """ get the corresponding mmidx path of the multireference dataset """
        if not os.path.isdir(CONFIG['apollo_reference_db_dir']):
            msg = "can't find the configured 'apollo_reference_db_dir': %s" % CONFIG['apollo_reference_db_dir']
            self.parser.error(msg)
        # validate multireference dataset at the subdirectory level
        try:
            validate_reference_dataset(CONFIG['apollo_reference_db_dir'])
        except FileNotFoundError as e:
            msg = "can't find required subdirectory in 'apollo_reference_db_dir': %s" % str(e)
            self.parser.error(''+msg)

        # validate multireference dataset at the actual file content level
        dfasm = read_reference_assembly_df()
        try:
            validate_reference_dataset(CONFIG['apollo_reference_db_dir'], dfasm)
        except FileNotFoundError as e:
            msg = "can't find required file in 'apollo_reference_db_dir': %s" % str(e)
            self.parser.error(''+msg)

        # now we can define the (existing) identity species mmidx
        mmidx = get_identify_species_mmidx_relpath(dfasm)
        mmidx = os.path.join(CONFIG['apollo_reference_db_dir'], "mmidx", mmidx)
        return mmidx

    def validate_apollo_reference(self,args=argparse.Namespace) -> bool:
        """ TODO validate ... """
        if args.custom_reference_fasta:
            if args.skip_reference_selection:
                self.identify_species_index = None
            else:
                # although external fasta provided, request is there to match-ref.py for informative purposes
                self.identify_species_index = self._get_mmidx_path()

            # TODO: custom reference needs to get copied at some point in time!
            self.species_reference = os.path.abspath(args.custom_reference_fasta)
            self.forced_species = True

        elif not args.custom_reference_tsv:
            # vanilla: run with standard multireference dataset
            self.identify_species_index = self._get_mmidx_path()

        elif args.custom_reference_tsv:
            explanation = """
                The option --custom-reference-tsv <TSV> is not supported (yet), but is here for explanatory purposes.
    
                When provided, the TSV should be validated that it is readable by:
                    df = read_reference_species_df()
                    ProvidedSchema.validate(df, lazy=True)
    
                This flag should (eventually ...) trigger:
                    - apollo-reference to kick in
                    - downloading all required genomes + mitochondria + taxonomical info (from NCBI)
                    - build reference fasta's and build species identification index
                    - deposit its output "somewhere"
                    - oh, and by the way, of course all error free ;-)
    
                This is for now a far-future scenario that is not being implemented.
                Use-case would be a (custom) clade group or a related group of species one would want to use.
                For now, this can be solved by:
                    - running the apollo-reference pipeline manually
                    - saving the complete output directory to a place accessible where compute will happen
                    - and provide this directory to --custom-reference-dataset </PATH/TO/REFERENCE_DIRECTORY>
    
                """
            print(explanation)
            sys.exit(0)


    def do_combinatorial_argument_validation(self,args=argparse.Namespace) -> bool:
        """ Combinatorial argument validation post parse_args() which can't be solved directly """

        # list of any of the --forced-*** flags for configured references
        reference_flags = [
            args.forced_species,
            args.forced_accession,
            args.forced_clade,
            args.forced_cladegroup,
            args.custom_reference_fasta
        ]


        # Namespace attributes can't be used together;
        # most of these are typically handled via parser.add_mutually_exclusive_group()
        argparse_mutual_exclusive_pairs = [
            ( 'forced_clade', 'skip_multiclade_strain_mapping' ),
            ( 'custom_clade_fasta', 'skip_multiclade_strain_mapping' ),
            ( 'custom_reference_fasta', 'skip_multiclade_strain_mapping' ),
        ]

        # Namespace attributes both required when being used
        argparse_mutual_required_pairs = []

        # second Namespace attribute requires first attribute to be set
        argparse_conditionally_required_pairs = [
            ( 'custom_reference_fasta', 'custom_reference_mitochondrion_fasta' ),
            ( 'custom_reference_fasta', 'custom_reference_mitochondrion_accessions' ),
            ( 'custom_clade_fasta', 'custom_clade_mitochondrion_fasta' ),
            ( 'custom_clade_fasta', 'custom_clade_mitochondrion_accessions' ),
        ]

        for _first, _second in argparse_mutual_exclusive_pairs:
            if getattr(args,_first) and getattr(args,_second):
                _first_option = get_flags_from_dest(self.parser,_first)[0]
                _second_option = get_flags_from_dest(self.parser,_second)[0]
                msg = "mutually exclusive: can't use %s in combination with %s" % (_first_option,_second_option)
                self.parser.error(msg)

        for _first, _second in argparse_mutual_required_pairs:
            if ( getattr(args,_first) and getattr(args,_second) ) != ( None, None ):
                _first_option = get_flags_from_dest(self.parser,_first)[0]
                _second_option = get_flags_from_dest(self.parser,_second)[0]
                msg = "mutually required: %s and %s needs to be specified together" % (_first_option,_second_option)
                self.parser.error(msg)

        for _first, _second in argparse_conditionally_required_pairs:
            if not getattr(args,_first) and getattr(args,_second):
                _first_option = get_flags_from_dest(self.parser,_first)[0]
                _second_option = get_flags_from_dest(self.parser,_second)[0]
                msg = "conditionally required: %s needs %s" % (_second_option,_first_option)
                self.parser.error(msg)

        # re'validate' default value for kraken.db_dir
        args.db_dir = validate_kraken_db(args.db_dir)

        if args.db_dir == None and not args.skip_kraken:
            msg = "can't find default kraken db path; consider using --skip-kraken [%s]" % CONFIG['kraken_db']
            self.parser.error(msg)

        elif args.skip_reference_selection and list(set(reference_flags))  == [None]:
            msg = "--skip_reference_selection requires one of -s/-a/-c/--clg/--custom-reference"
            self.parser.error(msg)

        elif args.custom_reference_fasta and args.custom_reference_mitochondrion_accessions:
            # TODO validate that args.custom_reference_mitochondrion_accessions occur in this fasta
            status = are_accessions_present_in_fasta(args.custom_reference_fasta,args.custom_reference_mitochondrion_accessions,warn=True)
            if not status:
                msg = "some --custom-reference-mitochondrion-accessions not found in fasta"
                self.parser.error(msg)

        elif args.custom_reference_fasta and args.custom_reference_mitochondrion_fasta:
            # TODO validate that accession(s) in args.custom_reference_mitochondrion_fasta don't occur in this fasta
            if args.custom_reference_fasta == args.custom_reference_mitochondrion_fasta:
                msg = "--custom-reference and --custom-reference-mitochondrion-fasta are the same!"
                self.parser.error(msg)
            status = are_no_accession_overlap_in_fasta(args.custom_reference_fasta,args.custom_reference_mitochondrion_fasta,warn=True)
            if not status:
                msg = "accession overlap in --custom-reference and --custom-reference-mitochondrion-accessions"
                self.parser.error(msg)

        elif args.extra_reference_fasta:
            # TODO validate that accession(s) in extra_reference_fasta don't occur in vanilla reference-tsv database
            # this is rather tricky: at this point in time, we've no access yet to the species selection database
            raise NotImplementedError("args.extra_reference_fasta")

        # in the scenario of --ISO(-species) flags, disallow most other flags to be defined
        elif args.ISO or [ k for k,v in args.__dict__.items() if k.startswith("ISO_") and v == True ] != []:
            assigned_reference_args = [ k for k, v in args.__dict__.items() if k.find("_reference_") > 0 and v not in (False,None) ]
            if assigned_reference_args != []:
                msg = f"--ISO and --ISO<-species> don't allow other pipeline-modifying arguments to be set {assigned_reference_args}"
                self.parser.error(msg)
            args.whatever = 100
            if args.ISO_Cauris:
                # hard-coded parameterization of analyses of C.auris clade samples
                pass

        # all combinations validated successfully!
        return True

    def do_forced_s_a_c_clg_argument_transformation(self,args=argparse.Namespace) -> bool:
        """ Transform --(forced-)accession, clade or cladegroup to the corresponding Namespace attribures

            - sets species_reference (fasta) and forced_species (just its name)
            - optionally sets clade_reference (fasta) forced_clade (just its name)

        """

        if ( args.forced_species or args.forced_accession or args.forced_clade or args.forced_cladegroup ):
            # !important! not read_reference_species_df() since we need reference assembly paths
            df = read_reference_species_and_assembly_df()
            if args.forced_species:
                # Link back species name corresponding fasta;
                # Mind that species name can occur multiple times in case of multi-clade references
                # With the example of candida_auris, default is to refer to the "central" is_primary=True reference
                index = df[((df.is_primary == True) | (df.is_primary.isnull()))][['species','fasta']].to_records(index=False).tolist()
                fasta = [ fasta for (species,fasta) in index if ( slugify_str_choice_type(species) == args.forced_species ) ][0]
                self.species_reference = os.path.join(CONFIG['apollo_reference_db_dir'],"refs",fasta)
                self.forced_species = args.forced_species
                return True
            elif (args.forced_accession or args.forced_clade or args.forced_cladegroup):
                if args.forced_accession:
                    _selected = df[df.reference==args.forced_accession]
                elif args.forced_cladegroup:
                    _selected = df[ ( (df.cladegroup==args.forced_cladegroup) & (df.is_primary == True) ) ]
                elif args.forced_clade:
                    # forced_clade sets a forced clade_reference AND its corresponding species_reference
                    _selected = df[df.clade==args.forced_clade.replace("_"," ")]
                    row = next(_selected.itertuples(index=False, name='Row'))
                    self.clade_reference = os.path.join(CONFIG['apollo_reference_db_dir'],"refs",row.fasta)
                    self.forced_clade = args.forced_clade
                    # now switch to the primary cladegroup's accession/reference
                    _selected = df[((df.cladegroup==row.cladegroup) & (df.is_primary == True) ) ]

                # in all these cases, redirect to the corresponding forced_species
                row = next(_selected.itertuples(index=False, name='Row'))
                args.forced_species = slugify_str_choice_type(row.species)
                # assign self.species_reference in any forced_**** scenario
                self.species_reference = os.path.join(CONFIG['apollo_reference_db_dir'],"refs",row.fasta)
                self.forced_species = args.forced_species
                return True
        else:
            return False

    def do_exterior_fasta_argument_transformation(self,args=argparse.Namespace) -> bool:
        """ redefine some arguments in the scenario of externally provided reference data """
        import warnings
        warnings.warn("do_exterior_fasta_argument_transformation need implementation!")
        return True

    def _parse_args(self) -> argparse.Namespace:
        args = super()._parse_args()

        # extra dataclass attributes of ApolloPipeline
        self.identify_species_index: Optional[Path] = None
        self.forced_species: Optional[str] = None
        self.forced_clade: Optional[str] = None
        self.species_reference: Optional[Path] = None
        self.clade_reference: Optional[Path] = None
        self.skip_reference_selection: Optional[bool] = args.__dict__['skip_reference_selection']
        self.skip_multiclade_strain_mapping: Optional[bool] = args.__dict__['skip_multiclade_strain_mapping']
        self.exterior_fasta: Optional[bool] = args.__dict__['custom_reference_fasta'] != None

        # TODO: handle mitochondrial accessions "somewhere"
        # - we can read from the "validate_apollo_reference" dataframe --> dict config
        # - we can decide, when copying multidb genome, to copy the mitochondrion.fal file alongside
        #   - but there are genomes without, so the "remove mt variant calls" has tiny extra complication
        #   - and there's the "external genome" scenario ...

        # do combinatorial checks and set reference fasta's to --forced-***
        self.validate_apollo_reference(args)
        self.do_combinatorial_argument_validation(args)
        self.do_forced_s_a_c_clg_argument_transformation(args)
        self.do_exterior_fasta_argument_transformation(args)

        self.kraken_db_dir: Path = args.db_dir
        # tool-specific parameter configuration
        self.mean_quality_threshold: int = args.__dict__['ithinkmultiqc.icannotfindthemptflag.mean_quality_threshold']
        self.window_size: int = args.__dict__['map_clean_reads.picard.window_size']
        self.min_read_length: int = args.__dict__['map_clean_reads.picard.minimum_length']

        # pipeline configuration
        self.time_limit: int = args.time_limit


        return args


    def setup(self) -> None:
        super().setup()

        if self.snakemake_args["use_singularity"]:
            self.snakemake_args["singularity_args"] = " ".join(
                [
                    self.snakemake_args["singularity_args"],
                    f"--bind {self.kraken_db_dir}:{self.kraken_db_dir}",
                ]  # paths that singularity should be able to read from can be bound by adding to the above list
            )

        # Change default time_limit to 180, or keep time_limit from the command line if > 180
        if self.time_limit < 300:
            self.time_limit = 300

        # DEPRECATED attributes (old Apollo pipeline with hard-coded references cauris and aspfum
        # select a reference based on species:
        # self.ref_dir = "/mnt/db/apollo/mapping/"
        # self.species == "candida_auris":
        # self.reference = Path("/mnt/db/apollo/mapping/candida_auris/GCA_002759435_3_genomic.fna")
        # self.custom_reference = /path/tp/fasta

        if self.species_reference is not None:
            print(f"# Running pipeline for species '{self.forced_species}' with reference: {self.species_reference}.")
        if self.clade_reference is not None:
            print(f"# Running pipeline for species '{self.forced_species}' with clade: {self.clade_reference}.")

        with open(
            Path(__file__).parent.joinpath("config/pipeline_parameters.yaml")
        ) as f:
            parameters_dict = yaml.safe_load(f)
        self.snakemake_config.update(parameters_dict)

        self.user_parameters = {
            "input_dir": str(self.input_dir),
            "output_dir": str(self.output_dir),
            "exclusion_file": str(self.exclusion_file),
            "mean_quality_threshold": int(self.mean_quality_threshold),
            "window_size": int(self.window_size),
            "min_read_length": int(self.min_read_length),
            "use_singularity": str(self.snakemake_args["use_singularity"]),
            "time-limit": str(self.time_limit),
            # selected or provided custom references (paths and names)
            "species_reference": str(self.species_reference),
            "clade_reference": str(self.clade_reference),
            "forced_species": str(self.forced_species),
            "forced_clade": str(self.forced_clade),
            "exterior_fasta": str(self.exterior_fasta),
            # database paths
            "reference_genomes_dir": str(os.path.join(CONFIG['apollo_reference_db_dir'],"refs")),
            "identify_species_index": str(self.identify_species_index),
            "kraken_db_dir": str(self.kraken_db_dir),
            # main switches that bypass/skip some parts of the snakemake workflow
            "skip_reference_selection": str(self.skip_reference_selection),
            "skip_multiclade_strain_mapping": str(self.skip_multiclade_strain_mapping),
        }


if __name__ == "__main__":
    main()

