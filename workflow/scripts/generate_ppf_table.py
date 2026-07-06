#!/usr/bin/env python3
from __future__ import annotations

""" Generate a PPF table for a specific SciPy statistical model and its parameters."""

# Python Imports
import numpy as np
import scipy.stats as stats
import pandas as pd
import yaml
from typing import Union, List, TYPE_CHECKING
import argparse

if TYPE_CHECKING:
    import numpy as np

ths = [0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.10, 0.25]
ths += [(1.0 - th) for th in reversed(ths)]
PPF_TABLE_DEFAULTS = ths

def generate_ppf_table(model_name: str, params: tuple, p_values: Union[np.ndarray|List[float]|None]=None, num_digits:Union[int|None]=None) -> pd.DataFrame:
    __doc__
    try:
        model = getattr(stats, model_name)
    except AttributeError:
        print(f"Error: Model '{model_name}' is no scipy.stats supported model.")
        return

    p_values = p_values or PPF_TABLE_DEFAULTS

    def _round(value:float,num_digits=num_digits) -> float:
        if num_digits is None:
            return value
        elif num_digits == 0:
            return int(value)
        else:
            return round(value, num_digits)

    # Calculate PPF values given each chance
    data = [ (p,_round(model.ppf(p, *params))) for p in p_values ]
    df = pd.DataFrame(data,columns=('p','threshold'))
    return df


def parse_p_values(p_str: str) -> list[float]:
    """Convert a comma-separated string of probabilities into a list of floats."""
    try:
        # Split by comma and convert each element to a float
        p_list = [float(p.strip()) for p in p_str.split(",")]

        # Validate that all probabilities are strictly between 0 and 1
        if any(p <= 0 or p >= 1 for p in p_list):
            raise argparse.ArgumentTypeError(
                "All p-values must be strictly between 0 and 1 (exclusive)."
            )
        return p_list
    except ValueError:
        raise argparse.ArgumentTypeError(
            "p-values must be a comma-separated list of numerical values."
        )


def get_parser() -> argparse.ArgumentParser:
    """ generate parser for this script """
    parser = argparse.ArgumentParser(
        description="Calculate a PPF (Percent Point Function) table for a given SciPy distribution model."
    )

    # Positional arguments
    parser.add_argument(
        "model_name",
        type=str,
        help="Name of the SciPy distribution model (e.g., 'norm', 't', 'weibull_min').",
    )
    parser.add_argument(
        "params",
        type=float,
        nargs="+",
        help="Space-separated parameters required for the chosen model (e.g., mean and std for 'norm').",
    )

    # Optional arguments
    parser.add_argument(
        "--p_values",
        type=parse_p_values,
        default=PPF_TABLE_DEFAULTS,
        help="Comma-separated list of explicit p-values (e.g., '0.05,0.5,0.95'). Defaults to a broad standard range.",
    )

    SUPPORTED_FORMATS = ["tsv", "csv", "yaml", "dict"]
    parser.add_argument(
        "--outstyle",
        type=str,
        choices=SUPPORTED_FORMATS,
        default="tsv",
        help="Output format style; choose from %s" % ", ".join(SUPPORTED_FORMATS),
    )

    parser.add_argument(
        "--round",
        type=int,
        default=None,
        choices=[0, 1, 2, 3, 4, 5, None],
        help="Round calulated PPF thresholds to at most N digits (default: not)",
    )


    return parser

def main() -> None:
    """ """
    parser = get_parser()
    args = parser.parse_args()

    # Run the generator with parsed arguments
    df = generate_ppf_table(
        model_name=args.model_name,
        params=args.params,
        p_values=args.p_values,
        num_digits=args.round,
    )

    if args.outstyle == "yaml":
        d = dict(df.itertuples(index=False, name=None))
        print(yaml.dump({"ppf_table": d }, sort_keys=False).rstrip("\n"))
    elif args.outstyle == "dict":
        print(dict(df.itertuples(index=False, name=None)))
    elif args.outstyle in ("tsv","csv"):
        sep = {"tsv":"\t","csv":","}.get(args.outstyle)
        print(df.to_csv(sep=sep,index=False).rstrip("\n"))
    else:
        raise ValueError("outstyle: '%s' is not supported." % args.outstyle)

if __name__ == "__main__":
    main()

