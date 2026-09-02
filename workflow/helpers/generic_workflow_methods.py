""" Placeholder module for generic workflow methods (see rivm-ids-swc-snakeutils).

For now, this code is kept here as pseudocode. It's currently not being used in apollo-mapping.

"""

# Imports
import inspect
import re
import sys
from typing import Literal, List, Tuple, Union, Dict, Any
from snakemake.workflow import Workflow

def get_rule_name() -> str:
    """
    Return the name of the closest Snakemake rule above the current line.

    This function inspects the Python call stack to access the workflow
    source code stored in the calling frame and determines which rule
    decorator appears closest above the current execution line.

    It searches for rule decorators of the form::

        @workflow.rule(name="rule_name", lineno=N)

    and returns the rule whose ``lineno`` is the largest value less than
    or equal to the current execution line.

    Returns
    -------
    str
        The name of the closest matching rule. If the rule name cannot
        be inferred, the fallback value ``"unknown_rule"`` is returned.


    Notes
    -----
    The function inspects two frames up the call stack to access the rule definition context.
    The frame reference is deleted afterward to avoid reference cycles when using the `inspect` module.
    """
    fallback = "unknown_rule"
    frame = inspect.currentframe()
    try:
        caller = frame.f_back if frame else None
        rule_frame = caller.f_back if caller else None
        code = getattr(rule_frame.f_locals, "get", lambda k, d=None: d)("code", "") if rule_frame else ""
        lineno = getattr(caller, "f_lineno", None) if caller else None

        if not code or not isinstance(lineno, int):
            return fallback

        # Name of the closest matching rule
        return max(
            (
                (m.group(1), int(m.group(2)))
                for m in re.finditer(r"@workflow\.rule\(name=['\"]([^'\"]+)['\"],\s*lineno=(\d+)", code)
                if int(m.group(2)) <= lineno
            ),
            default=(fallback, 0),
            key=lambda x: x[1],
        )[0]

    finally:
        del frame


def check_rule_environments(
        workflow: Workflow,
        env_type: Literal["conda", "container"] = "conda",
        strict: bool = False
) -> List[Tuple[str, str]]:
    """
    Inspects all rules in the provided Snakemake workflow to ensure they have defined
    a valid, non-empty conda environment or container image.

    Skipped rules include the 'all' rule, all localrules, any rule with 'base' in the
    underscore-splitted name, or rules opting out via params (e.g., 'threads_not_needed').

    Args:
        workflow: The global Snakemake workflow instance object.
        env_type: The environment type to check ("conda" or "container").
        strict: If True, raises a ValueError if any rule is missing the directive.
                If False, prints warnings to stderr and returns the missing entries.

    Returns:
        A list of tuples containing (rule_name, snakefile_path) for all rules
        missing the specified environment directive.
    """

    local_rules = getattr(workflow, "_localrules", set())
    missing_envs: List[Tuple[str, str]] = []

    for rule in workflow.rules:
        # Skip standard management rules, localrules or rules matching the 'base' naming convention
        if rule.name == "all" or rule.name in local_rules:
            continue
        elif "base" in rule.name.split("_"):
            continue
        elif hasattr(rule, "params"):
            # Safely check if 'conda_not_needed' or 'container_not_needed' is truthy
            bypass_flag = f"{env_type}_not_needed"
            if getattr(rule.params, bypass_flag, False) is True:
                continue

        is_missing = False

        if env_type == "conda":
            env_val = getattr(rule, "conda_env", None)
        elif env_type == "container":
            env_val = getattr(rule, "container_img", None)
        else:
            raise ValueError(f"Invalid env_type '{env_type}'. Choose 'conda' or 'container'.")

        if env_val is None or (isinstance(env_val, str) and not env_val.strip()):
            is_missing = True

        if is_missing:
            missing_envs.append((rule.name, rule.snakefile))

    if missing_envs:
        error_messages = [
            f"Rule '{r_name}' in file '{s_file}' is missing a valid {env_type} directive."
            for r_name, s_file in missing_envs
        ]

        if strict:
            raise ValueError(
                f"Strict environment check failed! The following rules lack a {env_type} environment:\n"
                + "\n".join(error_messages)
            )
        else:
            for msg in error_messages:
                print(f"WARNING: {msg}", file=sys.stderr)

    return missing_envs

def check_rule_directives(
        workflow: Workflow,
        directives: List[str],
        strict: bool = False
) -> List[Tuple[str, str, str]]:
    """
    Inspects all rules in the provided Snakemake workflow to ensure they have defined
    the specified directives or sub-properties (e.g., 'resources.mem_gb').

    Skipped rules include the 'all' rule, any rule with 'base' in the
    underscore-splitted name, or rules opting out via params (e.g., 'threads_not_needed').

    Example usage in your Snakefile (rather at the foot it):

        check_rule_directives(
            workflow=workflow,
            directives=["message", "threads", "resources.mem_gb", "log"],
            strict=False
        )

    Args:
        workflow: The global Snakemake workflow instance object.
        directives: A list of directive strings to check (e.g., ['threads', 'resources.mem_gb']).
        strict: If True, raises a ValueError if any rule is missing a directive.
                If False, prints warnings to stderr and returns the missing entries.

    Returns:
        A list of tuples containing (rule_name, snakefile_path, missing_directive)
        for all rules missing any of the specified directives.
    """
    missing_entries: List[Tuple[str, str, str]] = []

    for rule in workflow.rules:
        # Skip standard management rules or rules matching the 'base' naming convention
        if rule.name == "all":
            continue
        if "base" in rule.name.split("_"):
            continue

        for directive in directives:
            # Check if the rule explicitly opts out of this specific directive via params
            # e.g., if checking 'resources.mem_gb', it looks for 'resources_mem_gb_not_needed'
            if hasattr(rule, "params"):
                bypass_flag = f"{directive.replace('.', '_')}_not_needed"
                if getattr(rule.params, bypass_flag, False) is True:
                    continue

            is_missing = False
            parts = directive.split(".")
            base_directive = parts[0]

            # 1. Check the top-level directive attribute on the rule object
            attr_val = getattr(rule, base_directive, None)

            if attr_val is None:
                is_missing = True

            # 2. If it's a nested attribute (like resources.mem_gb)
            elif len(parts) > 1:
                sub_property = parts[1]

                # Snakemake resources/log/input usually act like dicts or custom objects
                if isinstance(attr_val, dict):
                    if sub_property not in attr_val or attr_val[sub_property] is None:
                        is_missing = True
                elif hasattr(attr_val, sub_property):
                    sub_val = getattr(attr_val, sub_property)
                    if sub_val is None:
                        is_missing = True
                else:
                    is_missing = True

            # 3. If found but it's an empty string or empty structure
            else:
                if isinstance(attr_val, str) and not attr_val.strip():
                    is_missing = True
                elif isinstance(attr_val, (list, dict, set)) and not attr_val:
                    is_missing = True

            if is_missing:
                missing_entries.append((rule.name, rule.snakefile, directive))

    if missing_entries:
        error_messages = [
            f"Rule '{r_name}' in file '{s_file}' is missing the '{dir_name}' directive."
            for r_name, s_file, dir_name in missing_entries
        ]

        if strict:
            raise ValueError(
                "Strict directive check failed! The following rules lack required directives:\n"
                + "\n".join(error_messages)
            )
        else:
            for msg in error_messages:
                print(f"WARNING: {msg}", file=sys.stderr)

    return missing_entries

def assign_mem_gb_default(
    workflow: Workflow,
    default_mem_gb: int) -> None:
    """ add resources['mem_gb'] to each rule that lacks the key 'mem_gb' in the rule.resources directive """
    local_rules = getattr(workflow, "_localrules", set())

    for rule in workflow.rules:
        # Skip standard management rules, localrules or rules matching the 'base' naming convention
        if rule.name == "all" or rule.name in local_rules:
            continue
        elif "base" in rule.name.split("_"):
            continue
        elif hasattr(rule,"resources") and ( hasattr(rule.resources, "mem_gb") or "mem_gb" in getattr(rule,"resources") ):
            continue
        print("assign_mem_gb_default:",(rule.name,default_mem_gb))
        rule.resources['mem_gb'] = default_mem_gb
