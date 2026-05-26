from __future__ import annotations

import os
from typing import TYPE_CHECKING
import FastaValidator
import tinyfasta

if TYPE_CHECKING:
    import argparse
    from pathlib import Path

def validate_fasta_file(arg:Path) -> Path:
    """ Validate that the provided argument is a valid fasta file.

    Use as: parer.add_argument("fasta",type=as_argparse_type(validate_fasta_file))
    """
    # https://pypi.org/project/py-fasta-validator/
    accepted_fasta_validator_return_codes = [0,2,4]
    if not os.path.isfile(str(arg)):
        msg = "not an existing (fasta) file: {}".format(arg)
        raise IOError(msg)
    return_code = FastaValidator.fasta_validator(arg)
    if return_code not in accepted_fasta_validator_return_codes:
        msg = f"invalid fasta file [py-fasta-validator: code={return_code}] {str(arg)}"
        raise ValueError(msg)
    return arg

def are_accessions_present_in_fasta(fasta:Path,accessions:list,warn=False) -> bool:
    """ Validate that all the accessions provided are present in the fasta file. """
    fa = tinyfasta.FastaParser(fasta)
    DESCRIPTIONS = [ str(rec.description).lstrip(">") for rec in fa ]
    ACCESSIONS = [ descr.split()[0] for descr in DESCRIPTIONS ]
    if not ACCESSIONS or not accessions:
        # No accessions in fasta or no accessions asked for.
        # This is considered incorrect behaviour of this function
        return False
    # convert to set of accessions and discard any shared one
    accessions = set(accessions)
    accessions.difference_update(DESCRIPTIONS)
    accessions.difference_update(ACCESSIONS)
    accessions = set([ acc.split()[0] for acc in accessions ])
    accessions.difference_update(ACCESSIONS)
    if not accessions:
        return True
    elif warn:
        import warnings
        for acc in accessions:
            warnings.warn(f"accession not in {fasta}: '{acc}'")
    return False


def are_accessions_absent_in_fasta(fasta:Path,accessions:list,warn=False) -> bool:
    """ Validate that none the accessions provided are present in the fasta file. """
    fa = tinyfasta.FastaParser(fasta)
    DESCRIPTIONS = [ str(rec.description).lstrip(">") for rec in fa ]
    ACCESSIONS = [ descr.split()[0] for descr in DESCRIPTIONS ]
    if not ACCESSIONS or not accessions:
        # No accessions in fasta or no accessions asked for.
        # This is considered incorrect behaviour of this function
        return False
    # convert to set of accessions and look for any shared accession
    accessions = set(accessions)
    shared = set()
    for _iter in (0,1):
        for _case in (DESCRIPTIONS, ACCESSIONS):
            if accessions.intersection(_case):
                shared.update(accessions.intersection(_case))
                accessions.difference_update(_case)
        # prepare final iteration: convert accessions to whitespace-splitted counterparts;
        # accessions were basically provided as 'descriptions'
        accessions = set([ acc.split()[0] for acc in accessions ])
    if not shared:
        return True
    elif warn:
        import warnings
        for acc in shared:
            warnings.warn(f"accession in {fasta}: '{acc}'")
    return False

def are_no_accession_overlap_in_fasta(fasta:Path,other:Path,warn=False) -> bool:
    """ Validate that none the accessions in two provided fasta files have identical names. """
    fa = tinyfasta.FastaParser(other)
    accessions = [ str(rec.description).lstrip(">") for rec in fa ]
    return are_accessions_absent_in_fasta(fasta,accessions,warn=warn)
