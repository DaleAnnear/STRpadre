from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from pathlib import Path

import build_consensus


ROOT = Path(__file__).parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


def call(caller: str, one: str, two: str = ".", status: str = "valid_call") -> dict[str, str]:
    return {"caller": caller, "call_status": status, "allele1_length_bp": one, "allele2_length_bp": two, "allele1_copy_number": ".", "allele2_copy_number": ".", "read_support_allele1": "10", "read_support_allele2": "10", "source_vcf": caller}


def config() -> dict:
    return {"minimum_caller_support": 2, "absolute_allele_length_tolerance_bp": 5, "relative_allele_length_tolerance": 0.05, "minimum_read_support": 0, "filtered_calls_may_contribute": False, "haploid_diploid_policy": "require_same_ploidy", "consensus_statistic": "median"}


def test_minimum_cost_pairing_handles_reversed_unphased_alleles() -> None:
    pairs, cost = build_consensus.assign_alleles([33.0, 30.0], [30.0, 33.0])
    assert pairs == [(0, 1), (1, 0)]
    assert cost == 0
    consensus = build_consensus.consensus_for_rows([call("longtr", "30", "33"), call("atarva", "33", "30")], config())
    assert consensus["status"] == "concordant"
    assert sorted(next(value for value in consensus.values() if isinstance(value, list) and value and isinstance(value[0], float))) == [30.0, 33.0]


def test_motif_rotations_and_reverse_complements_are_equivalent() -> None:
    policy = "rotations_and_reverse_complements"
    assert build_consensus.canonical_motif("CAG", policy) == build_consensus.canonical_motif("AGC", policy)
    assert build_consensus.canonical_motif("CAG", policy) == build_consensus.canonical_motif("CTG", policy)


def test_haploid_and_highly_discordant_calls() -> None:
    haploid = build_consensus.consensus_for_rows([call("longtr", "30"), call("atarva", "31")], config())
    assert haploid["status"] == "concordant"
    discordant = build_consensus.consensus_for_rows([call("longtr", "30", "33"), call("atarva", "90", "93")], config())
    assert discordant["status"] == "discordant"


def test_end_to_end_consensus_retains_matrix_statuses(tmp_path: Path) -> None:
    sample_sheet = tmp_path / "samples.csv"
    sample_sheet.write_text("sample_id,platform,alignment,alignment_index\nS1,hifi,a.bam,a.bai\nS2,ont,b.cram,b.crai\n")
    norm = tmp_path / "normal.tsv.gz"
    columns = ["sample_id", "platform", "canonical_locus_id", "chr", "canonical_start", "canonical_end", "reference_build", "motif", "caller", "call_status", "allele1_length_bp", "allele2_length_bp", "allele1_copy_number", "allele2_copy_number", "read_support_allele1", "read_support_allele2", "source_vcf"]
    rows = []
    for caller, one, two in [("trgt", "30", "33"), ("longtr", "31", "33"), ("atarva", "30", "34")]:
        rows.append({"sample_id": "S1", "platform": "hifi", "canonical_locus_id": "L1", "chr": "chr1", "canonical_start": "100", "canonical_end": "130", "reference_build": "test-build", "motif": "CAG", "caller": caller, "call_status": "valid_call", "allele1_length_bp": one, "allele2_length_bp": two, "allele1_copy_number": "10", "allele2_copy_number": "11", "read_support_allele1": "9", "read_support_allele2": "9", "source_vcf": caller + ".vcf"})
    with gzip.open(norm, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    out = tmp_path / "out"
    subprocess.run([sys.executable, str(ROOT / "bin" / "build_consensus.py"), "--normalized-files", str(norm), "--locus-manifest", str(FIXTURES / "manifest.tsv"), "--samplesheet", str(sample_sheet), "--config", str(ROOT / "configs" / "consensus.yml"), "--callers", "trgt,longtr,atarva,strdust", "--output-dir", str(out)], check=True)
    with gzip.open(out / "caller_matrix.tsv.gz", "rt") as handle:
        matrix = list(csv.DictReader(handle, delimiter="\t"))
    assert any(row["sample_id"] == "S2" and row["caller"] == "trgt" and row["caller_status"] == "not_applicable_to_platform" for row in matrix)
    assert (out / "consensus.vcf").is_file()
