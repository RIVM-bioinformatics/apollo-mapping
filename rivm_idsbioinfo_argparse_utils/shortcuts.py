# Python Imports
from typing import Protocol, Any, List
import argparse

class ArgumentContainer(Protocol):
    """ Typehint needed in functions have a (derived) argparse.ArgumentParser as argument

    Supported are:
        - argparse.ArgumentParser
        - argparse._ArgumentGroup
        - argparse._MutuallyExclusiveGroup

    Use-cases of where and when to use this:
        - in (complex) argparse menu's
        - argparse with subparser, with shared arguments that can't be solved with parents=[...]
        - toolbox of re-useable argparse recipies
        - and in all cases: required when strict code typing is enabled

    Example how to use in your own code repo:

        def add_custom_argument(container: ArgumentContainer) -> None:
            container.add_argument(
                "--flag",
                action="store_true",
                help="your help text",
            )

    Since when offered an argparse._ArgumentGroup, the below code has an incorrect type hint

        def add_custom_argument(parser: parser.ArgumentParser) -> None:
            parser.add_argument(
                "--flag",
                action="store_true",
                help="your help text",
            )



    """
    def add_argument(self, *args: str, **kwargs: Any) -> argparse.Action:
        ...



def get_flags_from_dest(parser:argparse.ArgumentParser, dest_name:str) -> List[str]:
    """ Returns a list of command-line flags for a given Namespace attribute name."""
    for action in parser._actions:
        if action.dest == dest_name:
            return action.option_strings
    return []

def as_argparse_type(validator):
    """ convert any validator (with a single argument) into an argparse type conversion function

    Allows re-use (DRY) of validators as argparse type conversion function. Example usage:

    def validate_something(input:Any) -> Any:
        if input == None:
            raise ValueError("can be anything, but not None")
        # optionally, modify input (in)to a modified or generalized type
        return input

    parser.add_argument(
        "--my-argument",
        type=as_argparse_type(validate_something),
        default=None,
        help="....."
    )

    """
    def wrapper(value):
        try:
            return validator(value)
        except Exception as e:
            # Translate any Exception into an argparse error
            raise argparse.ArgumentTypeError(str(e))
    return wrapper


