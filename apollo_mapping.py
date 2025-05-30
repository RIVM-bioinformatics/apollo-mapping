"""
Apollo mapping, based on Juno template
Authors: Roxanne Wolthuis, Boas van der Putten
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium
            Surveillance (IDS), Bacteriologie (BPD)     
Date: 07-04-2023   
"""

from pathlib import Path
import pathlib
import yaml
import argparse
import sys
from dataclasses import dataclass, field
from juno_library import Pipeline
from typing import Optional
from version import __package_name__, __version__, __description__


def main() -> None:
    apollo_mapping = ApolloMapping()
    apollo_mapping.run()


@dataclass
class ApolloMapping(Pipeline):
    pipeline_name: str = __package_name__
    pipeline_version: str = __version__
    input_type: str = "fastq"
    species_options = ["candida_auris", "aspergillus_fumigatus"]

    def _add_args_to_parser(self) -> None:
        super()._add_args_to_parser()

        self.parser.description = (
            "Apollo mapping pipelines for reference mapping analysis of fungal genomes."
        )

        self.add_argument(
            "-s",
            "--species",
            type=str.lower,
            metavar="STR",
            help=f"Species to use, choose from: {self.species_options}",
            required=True,
            dest="species",
            choices=self.species_options,
        )
        self.add_argument(
            "--reference",
            type=Path,
            metavar="FILE",
            dest="custom_reference",
            help="Reference genome to use default is chosen based on species argument, defaults per species can be found in: /mnt/db/apollo/mapping/[species]",
            required=False,
        )
        self.add_argument(
            "--db-dir",
            type=Path,
            default="/mnt/db/apollo/kraken_db_apollo",
            metavar="DIR",
            help="Kraken2 database directory (should include fungi!).",
        )
        self.add_argument(
            "-mpt",
            "--mean-quality-threshold",
            type=int,
            metavar="INT",
            default=28,
            help="Phred score to be used as threshold for cleaning (filtering) fastq files.",
        )
        self.add_argument(
            "-ws",
            "--window-size",
            type=int,
            metavar="INT",
            default=5,
            help="Window size to use for cleaning (filtering) fastq files.",
        )
        self.add_argument(
            "-ml",
            "--minimum-length",
            type=int,
            metavar="INT",
            default=50,
            help="Minimum length for fastq reads to be kept after trimming.",
        )

    def _parse_args(self) -> argparse.Namespace:
        args = super()._parse_args()

        # Optional arguments are loaded into self here
        self.db_dir: Path = args.db_dir
        self.mean_quality_threshold: int = args.mean_quality_threshold
        self.window_size: int = args.window_size
        self.min_read_length: int = args.minimum_length
        self.reference: Optional[Path] = None
        self.custom_reference: Path = args.custom_reference
        self.time_limit: int = args.time_limit
        self.species: str = args.species

        return args

    def setup(self) -> None:
        super().setup()

        if self.snakemake_args["use_singularity"]:
            self.snakemake_args["singularity_args"] = " ".join(
                [
                    self.snakemake_args["singularity_args"],
                    f"--bind {self.db_dir}:{self.db_dir}",
                ]  # paths that singularity should be able to read from can be bound by adding to the above list
            )

        # Change default time_limit to 180, or keep time_limit from the command line if > 180
        if self.time_limit < 300:
            self.time_limit = 300

        # select a reference based on species:
        # self.ref_dir = "/mnt/db/apollo/mapping/"
        if self.species == "candida_auris":
            self.reference = Path(
                "/mnt/db/apollo/mapping/candida_auris/GCF_003013715.1.fasta"
            )
        elif self.species == "aspergillus_fumigatus":
            self.reference = Path(
                # "/mnt/db/apollo/mapping/aspergillus_fumigatus/CEA10.fasta"
                "/mnt/db/apollo/mapping/aspergillus_fumigatus/GCF_000002655_1.fna"
            )

        if self.custom_reference is not None:
            print(
                "A reference genome was specified by the user, which may not be the default reference genome for this species."
            )
            self.reference = self.custom_reference

        print(f"Running pipeline for {self.species} with reference: {self.reference}.")

        with open(
            Path(__file__).parent.joinpath("config/pipeline_parameters.yaml")
        ) as f:
            parameters_dict = yaml.safe_load(f)
        self.snakemake_config.update(parameters_dict)

        self.user_parameters = {
            "input_dir": str(self.input_dir),
            "output_dir": str(self.output_dir),
            "exclusion_file": str(self.exclusion_file),
            "db_dir": str(self.db_dir),
            "mean_quality_threshold": int(self.mean_quality_threshold),
            "window_size": int(self.window_size),
            "min_read_length": int(self.min_read_length),
            "reference": str(self.reference),
            "use_singularity": str(self.snakemake_args["use_singularity"]),
            "time-limit": str(self.time_limit),
            "species": str(self.species),
        }


if __name__ == "__main__":
    main()
