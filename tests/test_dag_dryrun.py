import io
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar

import pandas as pd
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = PROJECT_ROOT / "tests" / "fixtures" / "fastq" / "illumina_PE_empty"


@dataclass(frozen=True)
class DagRuleExpectation:
    """Describe whether a Snakemake rule should appear and how many jobs it has.

    A count of zero expects the rule to be absent. Positive counts specify the
    per-sample multiplier by default; set ``per_sample`` false for a fixed total.
    """

    rule: str = field(metadata={"description": "Name of the Snakemake rule."})
    count: int = field(
        metadata={
            "description": "Expected job count; zero means the rule must be absent."
        }
    )
    per_sample: bool = field(
        default=True,
        metadata={
            "description": (
                "Whether a positive count is a multiplier of fixture pairs "
                "(default) or a fixed total job count."
            )
        },
    )

    def expected_count(self, sample_pairs: int) -> int:
        """Return the expected job count for the fixture's sample pairs."""
        return self.count * sample_pairs if self.per_sample else self.count


@dataclass(frozen=True)
class DagDryrunCase:
    """Bundle CLI arguments with shared default and case-specific DAG expectations."""

    workflow_args: tuple[str, ...] = field(
        metadata={
            "description": "CLI arguments passed to apollo_mapping.py for this workflow case."
        }
    )
    additional_rules: tuple[DagRuleExpectation, ...] = ()
    # DagRuleExpectation rules that should be valid for any Snakemake DAG dry-run case
    rules: ClassVar[tuple[DagRuleExpectation, ...]] = (
        DagRuleExpectation("all", count=1, per_sample=False),
        # Consider to add more in a harmonized pipeline landscape:
        #DagRuleExpectation("multiqc", count=1, per_sample=False),
        #DagRuleExpectation("collectresults", count=1, per_sample=False),
    )


def parse_dag_table_from_dryrun(dryrun_output: str) -> pd.DataFrame:
    """Parse Snakemake's job statistics table and omit its aggregate total row."""
    job_stats = dryrun_output.split("Job stats:", maxsplit=1)
    if len(job_stats) != 2:
        raise ValueError("Dry-run output does not contain a Snakemake job stats table")

    table = job_stats[1].split("\n\n", maxsplit=1)[0].strip()
    dag = pd.read_fwf(io.StringIO(table), skiprows=[1])
    if not {"job", "count"}.issubset(dag.columns):
        raise ValueError("Could not parse the Snakemake job stats table")
    return dag.loc[dag["job"] != "total"].reset_index(drop=True)


def expected_sample_pair_count() -> int:
    """Count fixture sample pairs, raising if either read mate is missing."""
    read1_samples = {
        path.name.removesuffix("_1.fq.gz")
        for path in FIXTURES.glob("*_1.fq.gz")
    }
    read2_samples = {
        path.name.removesuffix("_2.fq.gz")
        for path in FIXTURES.glob("*_2.fq.gz")
    }
    if not read1_samples or read1_samples != read2_samples:
        raise AssertionError("Fixture FASTQs must contain matching R1/R2 sample pairs")
    return len(read1_samples)



# DagRuleExpectation that are valid in each parameterization of the DAG dry-run tests
GENERALLY_EXPECTED_DAG_RULES = (
    DagRuleExpectation("assign_reference", count=1),
    DagRuleExpectation("forced_ref", count=1),
    DagRuleExpectation("bwa_mem", count=1),
    DagRuleExpectation("mason_genome", count=1, per_sample=False),
)

# should be absent in all DAGs except for --trigger-multiclade-masking-workflow
GENERALLY_UNEXPECTED_DAG_RULES = (
    DagRuleExpectation("bam_to_softclipped_bed", count=0),
)
@pytest.mark.parametrize(
    "case",
    [
        pytest.param(
            DagDryrunCase(
                workflow_args=("--species", "candida_auris", "--skip-kraken"),
                additional_rules=(
                    *GENERALLY_EXPECTED_DAG_RULES,
                    *GENERALLY_UNEXPECTED_DAG_RULES,
                    # negative test: this rule should not exist in the DAG
                    DagRuleExpectation("non_existing_rule", count=0),
                ),
            ),
            id="species",
        ),
        pytest.param(
            DagDryrunCase(
                workflow_args=("--clg", "cauris", "--skip-kraken"),
                additional_rules=(
                    *GENERALLY_EXPECTED_DAG_RULES,
                    *GENERALLY_UNEXPECTED_DAG_RULES,
                    DagRuleExpectation("match_ref", count=1),
                ),
            ),
            id="cladegroup",
        ),
        pytest.param(
            DagDryrunCase(
                workflow_args=(
                    "--species",
                    "candida_auris",
                    "--skip-reference-selection",
                    "--skip-kraken",
                ),
                additional_rules=(
                    *GENERALLY_EXPECTED_DAG_RULES,
                    *GENERALLY_UNEXPECTED_DAG_RULES,
                    # rule should not be in the DAG due to --skip-reference-selection
                    DagRuleExpectation("match_ref", count=0),
                ),
            ),
            id="skip-reference-selection",
        ),
        pytest.param(
            DagDryrunCase(
                workflow_args=(
                    "--species",
                    "candida_auris",
                    "--skip-multiclade-strain-mapping",
                    "--skip-kraken",
                ),
                additional_rules=(
                    *GENERALLY_EXPECTED_DAG_RULES,
                    *GENERALLY_UNEXPECTED_DAG_RULES,
                ),
            ),
            id="skip-multiclade-strain-mapping",
        ),
        pytest.param(
            DagDryrunCase(
                workflow_args=(
                    "--species",
                    "candida_auris",
                    "--trigger-multiclade-masking-workflow",
                ),
                additional_rules=(
                    *GENERALLY_EXPECTED_DAG_RULES,
                    DagRuleExpectation("bam_to_softclipped_bed", count=1),
                ),
            ),
            id="multiclade-masking",
        ),
    ],
)
def test_snakemake_dag_dryrun(tmp_path: Path, case: DagDryrunCase) -> None:
    """Run each configured CLI dry-run and assert its expected Snakemake DAG."""
    for name in ("Snakefile", "data", "envs", "files", "juno_mapping", "workflow"):
        source = PROJECT_ROOT / name
        (tmp_path / name).symlink_to(source, target_is_directory=source.is_dir())
    shutil.copytree(PROJECT_ROOT / "config", tmp_path / "config")
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    command = [
        sys.executable,
        str(PROJECT_ROOT / "apollo_mapping.py"),
        "--input",
        str(FIXTURES),
        "--output",
        str(output_dir),
        "--workdir",
        str(tmp_path),
        *case.workflow_args,
        # all these below needed to trigger class Pipeline to execute locally
        "--local",
        "--no-containers",
        "--snakemake-args",
        "cores=1",
        "nodes=1",
        "--dryrun",
    ]

    result = subprocess.run(
        command,
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )
    output = result.stdout + result.stderr

    assert result.returncode == 0, output
    dag = parse_dag_table_from_dryrun(output)
    job_counts = dag.set_index("job")["count"]
    expected_pairs = expected_sample_pair_count()

    for expectation in (*case.rules, *case.additional_rules):
        expected_count = expectation.expected_count(expected_pairs)
        if expected_count == 0:
            assert expectation.rule not in job_counts, output
        else:
            assert expectation.rule in job_counts, output
            assert job_counts[expectation.rule] == expected_count, output
