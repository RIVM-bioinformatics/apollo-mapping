import inspect
import re

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
