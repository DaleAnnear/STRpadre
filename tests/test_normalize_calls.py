from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


@pytest.mark.parametrize("caller", ["trgt", "longtr", "atarva", "strdust"])
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
    assert row["parser_schema_version"] == "1.0.0"


def test_filtered_and_haploid_calls_do_not_invent_second_allele(tmp_path: Path) -> None:
    vcf = tmp_path / "haploid.vcf"
    vcf.write_text("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t101\t.\tCAG\tCAGCAG\t.\tLOW_SUPPORT\tEND=130\tGT:FRB:SUP\t1:30:4\n")
    output, versions = tmp_path / "out.tsv.gz", tmp_path / "versions.yml"
    subprocess.run([sys.executable, str(ROOT / "bin" / "normalize_calls.py"), "--caller", "strdust", "--sample-id", "S1", "--platform", "ont", "--vcf", str(vcf), "--locus-manifest", str(FIXTURES / "manifest.tsv"), "--source-vcf", "native.vcf", "--output", str(output), "--versions", str(versions)], check=True)
    with gzip.open(output, "rt") as handle: row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["call_status"] == "filtered_call"
    assert row["allele2_length_bp"] == "."
