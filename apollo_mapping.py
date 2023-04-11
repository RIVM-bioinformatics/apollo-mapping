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

    def _add_args_to_parser(self) -> None:
        super()._add_args_to_parser()

        self.parser.description = "Apollo mapping pipelines for reference mapping analysis of fungal genomes."
        
        self.add_argument(
            "--reference",
            type=Path,
            metavar="FILE",
            help="Reference genome to use",
            required=True
        )
        self.add_argument(
            "--db-dir",
            type=Path,
            default="",
            metavar="DIR",
            help="Kraken2 database directory (should include fungi!).",
        )
        self.add_argument(
            "-mpt",
            "--mean-quality-threshold",
            type = int,
            metavar = "INT",
            default = 28,
            help = "Phred score to be used as threshold for cleaning (filtering) fastq files."
        )
        self.add_argument(
            "-ws",
            "--window-size",
            type = int,
            metavar = "INT",
            default = 5,
            help = "Window size to use for cleaning (filtering) fastq files."
        )
        self.add_argument(
            "-ml",
            "--minimum-length",
            type = int,
            metavar = "INT",
            default = 50,
            help = "Minimum length for fastq reads to be kept after trimming."
        )
        
    def _parse_args(self) -> argparse.Namespace:
        args = super()._parse_args()

        # Optional arguments are loaded into self here
        self.db_dir: Path = args.db_dir
        self.mean_quality_threshold: int = args.mean_quality_threshold
        self.window_size: int = args.window_size
        self.min_read_length: int = args.minimum_length
        self.reference: Path = args.reference

        return args
    
    # # Extra class methods for this pipeline can be defined here
    # def example_class_method(self):
    #     print(f"example option is set to {self.example}")

    def setup(self) -> None:
        super().setup()

        if self.snakemake_args["use_singularity"]:
            self.snakemake_args["singularity_args"] = " ".join(
                [
                    self.snakemake_args["singularity_args"]
                ] # paths that singularity should be able to read from can be bound by adding to the above list
            )

        # # Extra class methods for this pipeline can be invoked here
        # if self.example:
        #     self.example_class_method()

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
            "use_singularity": str(self.snakemake_args['use_singularity']),
        }


if __name__ == "__main__":
    main()
