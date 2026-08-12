from __future__ import annotations

import gzip

from argparse import Namespace
from pathlib import Path

import run_caller


def longtr_command(options: dict) -> list[str]:
    config = {
        "threads": 1,
        "phased_reads_expected": False,
        "additional_args": [],
        "options": options,
    }
    args = Namespace(
        alignment=Path("reads.bam"),
        reference=Path("reference.fa"),
        catalog=Path("catalog.bed"),
        sample_id="SAMPLE",
    )
    command, _, _ = run_caller.build_command("longtr", config, args, Path("output"))
    return command


def test_longtr_uses_sample_as_library_by_default() -> None:
    assert "--lib-from-samp" in longtr_command({})


def test_longtr_can_use_read_group_library_tags() -> None:
    assert "--lib-from-samp" not in longtr_command({"use_lb_tags": True})

def test_longtr_missing_dflankindel_header_is_repaired(tmp_path: Path) -> None:
    raw_vcf = tmp_path / "calls.original.vcf.gz"
    raw_text = (
        "##fileformat=VCFv4.3\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\t.\tA\tC\t.\tPASS\tEND=101\tGT:DFLANKINDEL\t0/1:1,2\n"
    )
    with gzip.open(raw_vcf, "wt") as handle:
        handle.write(raw_text)

    repaired_vcf = run_caller.repair_longtr_vcf_header(raw_vcf)

    assert repaired_vcf != raw_vcf
    with gzip.open(repaired_vcf, "rt") as handle:
        repaired_text = handle.read()
    assert "##FORMAT=<ID=DFLANKINDEL,Number=.,Type=String," in repaired_text
    assert repaired_text.index("##FORMAT=<ID=DFLANKINDEL") < repaired_text.index("#CHROM")
