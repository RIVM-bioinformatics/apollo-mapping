""" Placeholder module for generic workflow methods (see rivm-ids-swc-snakeutils).

For now, this code is kept here as pseudocode. It's currently not being used in apollo-mapping.

"""

# Imports
import inspect
import re
import sys
from typing import Literal, List, Tuple
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


import sys
from typing import Literal, List, Tuple
from snakemake.workflow import Workflow


def check_rule_environments(
        workflow: Workflow,
        env_type: Literal["conda", "container"] = "conda",
        strict: bool = False
) -> List[Tuple[str, str]]:
    """
    Inspects all rules in the provided Snakemake workflow to ensure they have defined
    a valid, non-empty conda environment or container image, skipping rules explicitly
    marked as not needing one via the params directive,
    or any rule with 'base' in the underscore-splitted name.

    Args:
        workflow: The global Snakemake workflow instance object.
        env_type: The environment type to check ("conda" or "container").
        strict: If True, raises a ValueError if any rule is missing the directive.
                If False, prints warnings to stderr and returns the missing entries.

    Returns:
        A list of tuples containing (rule_name, snakefile_path) for all rules
        missing the specified environment directive.
    """
    missing_envs: List[Tuple[str, str]] = []

    for rule in workflow.rules:
        # Check if the rule explicitly opts out using its name or given params
        # rule.params is an object that behaves like a dictionary/namespace
        if rule.name == "all":
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
