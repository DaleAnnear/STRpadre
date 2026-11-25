from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


@pytest.mark.parametrize("caller", ["trgt", "longtr", "atarva", "strdust", "strkit"])
def test_each_documented_parser_emits_common_schema(tmp_path: Path, caller: str) -> None:
    output = tmp_path / "calls.tsv.gz"
    versions = tmp_path / "versions.yml"
    command = [sys.executable, str(ROOT / "bin" / "normalize_calls.py"), "--caller", caller, "--sample-id", "S1", "--platform", "hifi", "--vcf", str(FIXTURES / f"{caller}.vcf"), "--locus-manifest", str(FIXTURES / "manifest.tsv"), "--source-vcf", f"raw/{caller}/S1/calls.vcf.gz", "--output", str(output), "--versions", str(versions)]
    subprocess.run(command, check=True)
    with gzip.open(output, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["canonical_locus_id"] == "L1"
    assert row["call_status"] == "valid_call"
    assert row["allele1_length_bp"] == "30"
    assert row["allele2_length_bp"] == "33"
    if caller == "strkit":
        assert (row["allele1_sequence"], row["allele2_sequence"]) == ("CAG" * 10, "CAG" * 11)
        assert (row["allele1_copy_number"], row["allele2_copy_number"]) == ("10", "11")
        assert (row["read_support_allele1"], row["read_support_allele2"]) == ("7", "8")
        assert row["phase_set"] == "42"
    assert row["parser_schema_version"] == "1.0.0"


@pytest.mark.parametrize(
    ("gt", "ploidy", "sd", "status", "lengths", "support"),
    [
        ("1|0", 2, "7,8", "valid_call", ("33", "30"), ("7", "8")),
        ("./1", 1, "7", "valid_call", ("33", "."), ("7", ".")),
        (".|1", 2, "7", "valid_no_call", (".", "33"), (".", "7")),
        ("1/1", 2, "15", "valid_call", ("33", "33"), (".", ".")),
        ("2/0", 2, "5,6", "valid_call", ("0", "30"), ("5", "6")),
    ],
)
def test_medaka_replacement_alleles(gt, ploidy, sd, status, lengths, support):
    from normalize_calls import VcfRecord, load_manifest, map_locus, parse_record

    record = VcfRecord(
        "chr1", 101, "chr1_100_130", "CAG" * 10, ["CAG" * 11, "<DEL>"], "PASS",
        {"rec": f"tr_chr1_100_130_pad_90_140_fwd_hap0_phased-set0_ploidy{ploidy}"},
        {"GT": gt, "SD": sd},
    )
    loci, _ = load_manifest(FIXTURES / "manifest.tsv")
    row = parse_record("medaka", record, map_locus("medaka", record, loci), "S1", "ont", "raw.vcf")
    assert row["call_status"] == status
    assert (row["allele1_length_bp"], row["allele2_length_bp"]) == lengths
    assert (row["read_support_allele1"], row["read_support_allele2"]) == support
    assert row["allele1_copy_number"] == row["allele2_copy_number"] == "."


def test_filtered_and_haploid_calls_do_not_invent_second_allele(tmp_path: Path) -> None:
    vcf = tmp_path / "haploid.vcf"
    vcf.write_text("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t101\t.\tCAG\tCAGCAG\t.\tLOW_SUPPORT\tEND=130\tGT:FRB:SUP\t1:30:4\n")
    output, versions = tmp_path / "out.tsv.gz", tmp_path / "versions.yml"
    subprocess.run([sys.executable, str(ROOT / "bin" / "normalize_calls.py"), "--caller", "strdust", "--sample-id", "S1", "--platform", "ont", "--vcf", str(vcf), "--locus-manifest", str(FIXTURES / "manifest.tsv"), "--source-vcf", "native.vcf", "--output", str(output), "--versions", str(versions)], check=True)
    with gzip.open(output, "rt") as handle: row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["call_status"] == "filtered_call"
    assert row["allele2_length_bp"] == "."
