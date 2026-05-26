""" copied from apollo-variant-typing/apollo_variant_typing.py / juno-variant-typing/juno_variant_typing.py """

# Python Imports
import argparse
from typing import Union, Callable

def check_number_within_range(
    minimum: float = 0, maximum: float = 1
) -> Union[Callable[[str], str], argparse.FileType]:
    """
    Creates a function to check whether a numeric value is within a range, inclusive.

    The generated function can be used by the `type` parameter in argparse.ArgumentParser.
    See https://stackoverflow.com/a/53446178.

    Args:
        value: the numeric value to check.
        minimum: minimum of allowed range, inclusive.
        maximum: maximum of allowed range, inclusive.

    Returns:
        A function which takes a single argument and checks this against the range.

    Raises:
        argparse.ArgumentTypeError: if the value is outside the range.
        ValueError: if the value cannot be converted to float.
    """

    def generated_func_check_range(value: str) -> str:
        value_f = float(value)
        if (value_f < minimum) or (value_f > maximum):
            raise argparse.ArgumentTypeError(
                f"Supplied value {value} is not within expected range {minimum} to {maximum}."
            )
        return str(value)

    return generated_func_check_range