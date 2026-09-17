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

def sex_aware_config() -> dict:
    return {
        "threads": 1,
        "phased_reads_expected": False,
        "additional_args": [],
        "options": {"haploid_chromosomes": []},
        "sex_chromosome_policy": {
            "enabled": True,
            "reference_build": "GRCh38",
            "chromosome_x": "chrX",
            "chromosome_y": "chrY",
            "pseudoautosomal_regions": [
                {"chromosome": "chrX", "start": 10000, "end": 2781479},
                {"chromosome": "chrX", "start": 155701382, "end": 156030895},
                {"chromosome": "chrY", "start": 10000, "end": 2781479},
                {"chromosome": "chrY", "start": 56887902, "end": 57217415},
            ],
        },
    }


def test_male_sex_chromosome_catalog_is_partitioned_by_par(tmp_path: Path) -> None:
    catalog = tmp_path / "longtr.bed"
    catalog.write_text(
        "chr1\t101\t130\tCAG\tAUTOSOME\n"
        "chrX\t154437173\t154437196\tGGC\tX_NON_PAR\n"
        "chrX\t155800001\t155800010\tGGC\tX_PAR2\n"
        "chrY\t2935976\t2935991\tCGG\tY_NON_PAR\n"
        "chrY\t56888001\t56888010\tCGG\tY_PAR2\n"
    )
    args = Namespace(catalog=catalog, sex="male", sample_id="S1")

    invocations, plan = run_caller.sex_aware_invocations("longtr", sex_aware_config(), args, tmp_path)

    assert plan == {"enabled": True, "sex": "male", "haploid_loci": 2, "diploid_loci": 3, "skipped_y_loci": 0}
    assert [invocation.label for invocation in invocations] == ["diploid", "haploid"]
    diploid, haploid = invocations
    assert "X_PAR2" in diploid.catalog.read_text() and "Y_PAR2" in diploid.catalog.read_text()
    assert "X_NON_PAR" not in diploid.catalog.read_text() and "Y_NON_PAR" not in diploid.catalog.read_text()
    assert "X_NON_PAR" in haploid.catalog.read_text() and "Y_NON_PAR" in haploid.catalog.read_text()
    assert haploid.config["options"]["haploid_chromosomes"] == ["chrX", "chrY"]
    assert diploid.config["options"]["haploid_chromosomes"] == []


def test_female_sex_chromosome_catalog_skips_y_loci(tmp_path: Path) -> None:
    catalog = tmp_path / "strdust.bed"
    catalog.write_text("chrX\t154437172\t154437196\nchrY\t2935975\t2935991\n")
    args = Namespace(catalog=catalog, sex="XX", sample_id="S1")

    invocations, plan = run_caller.sex_aware_invocations("strdust", sex_aware_config(), args, tmp_path)

    assert plan == {"enabled": True, "sex": "female", "haploid_loci": 0, "diploid_loci": 1, "skipped_y_loci": 1}
    assert [invocation.label for invocation in invocations] == ["diploid"]
    assert "chrX" in invocations[0].catalog.read_text()
    assert "chrY" not in invocations[0].catalog.read_text()


def test_sex_aware_catalog_requires_sex_metadata(tmp_path: Path) -> None:
    catalog = tmp_path / "longtr.bed"
    catalog.write_text("chrX\t154437173\t154437196\tGGC\tX_NON_PAR\n")
    args = Namespace(catalog=catalog, sex="", sample_id="S1")

    try:
        run_caller.sex_aware_invocations("longtr", sex_aware_config(), args, tmp_path)
    except run_caller.CallerError as error:
        assert "sex must be" in str(error)
    else:
        raise AssertionError("missing sex metadata accepted")

def test_sex_aware_vcfs_are_merged_then_coordinate_sorted(tmp_path: Path, monkeypatch) -> None:
    diploid = tmp_path / "calls.diploid.original.vcf"
    haploid = tmp_path / "calls.haploid.original.vcf"
    for path in (diploid, haploid):
        path.write_text("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
    commands: list[list[str]] = []
    monkeypatch.setattr(run_caller, "executable", lambda name: name)
    monkeypatch.setattr(run_caller, "run", lambda command, log, stdout=None: commands.append(command))

    run_caller.normalise_vcfs([diploid, haploid], tmp_path / "calls.vcf.gz", tmp_path / "caller.stderr.log", "strdust")

    assert [command[1] for command in commands if command[0] == "bcftools"] == ["sort", "sort", "concat", "sort"]
    assert commands[-1] == ["tabix", "-f", "-p", "vcf", str(tmp_path / "calls.vcf.gz")]
