from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import gzip
import yaml

import validate_inputs


def caller_config(caller: str, build: str = "test-build") -> dict:
    versions = {"trgt": "5.1.0", "longtr": "1.2", "atarva": "0.7.1", "strdust": "0.20.0"}
    platforms = {"trgt": ["hifi"], "longtr": ["hifi", "ont"], "atarva": ["hifi", "ont"], "strdust": ["hifi", "ont"]}
    return {"schema_version": "1.0", "caller": caller, "container": f"local/{caller}:{versions[caller]}", "tool_version": versions[caller], "reference_build": build, "platforms": platforms[caller], "native_catalog": "ignored", "threads": 1 if caller == "longtr" else 2, "resources": {"memory": "4 GB", "time": "1h"}, "phased_reads_expected": False, "additional_args": [], "options": {}}


def invoke(monkeypatch, args: list[str]) -> int:
    monkeypatch.setattr(sys, "argv", ["validate_inputs.py", *args])
    return validate_inputs.main()


def test_preflight_restricts_trgt_to_hifi(tmp_path: Path, monkeypatch) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\n" + "N" * 200 + "\n")
    (tmp_path / "ref.fa.fai").write_text("chr1\t200\t6\t200\t201\n")
    manifest = tmp_path / "loci.tsv"
    manifest.write_text((Path(__file__).parent / "fixtures" / "manifest.tsv").read_text())
    bam, bai, cram, crai = (tmp_path / "a.bam", tmp_path / "a.bam.bai", tmp_path / "b.cram", tmp_path / "b.cram.crai")
    for path in (bam, bai, cram, crai): path.touch()
    samples = tmp_path / "samples.csv"
    samples.write_text("sample_id,platform,alignment,alignment_index\nS_HIFI,hifi," + str(bam) + "," + str(bai) + "\nS_ONT,ont," + str(cram) + "," + str(crai) + "\n")
    catalogs = {}
    for caller, content in {"trgt": "chr1\t100\t130\tID=L1;MOTIFS=CAG;STRUC=<TR>\n", "longtr": "chr1\t101\t130\tCAG\tL1\n", "strdust": "chr1\t100\t130\n"}.items():
        catalogs[caller] = tmp_path / f"{caller}.bed"; catalogs[caller].write_text(content)
    atarva_plain = tmp_path / "atarva.bed"
    atarva_plain.write_text("chr1\t100\t130\tCAG\t3\n")
    atarva = tmp_path / "atarva.bed.gz"
    with gzip.open(atarva, chr(119) + chr(116)) as handle: handle.write(atarva_plain.read_text()); Path(str(atarva) + chr(46) + chr(116) + chr(98) + chr(105)).touch()
    catalogs["atarva"] = atarva
    configs = {}
    for caller in ("trgt", "longtr", "atarva", "strdust"):
        configs[caller] = tmp_path / f"{caller}.yml"; configs[caller].write_text(yaml.safe_dump(caller_config(caller)))
    consensus = tmp_path / "consensus.yml"
    consensus.write_text((Path(__file__).parents[1] / "configs" / "consensus.yml").read_text())
    out = tmp_path / "out"
    args = ["--samplesheet", str(samples), "--reference", str(reference), "--reference-fai", str(reference) + ".fai", "--locus-manifest", str(manifest), "--callers", "trgt,longtr,atarva,strdust", "--consensus-config", str(consensus), "--output-dir", str(out), "--skip-alignment-header-check"]
    for caller in ("trgt", "longtr", "atarva", "strdust"):
        args += [f"--{caller}-config", str(configs[caller]), f"--{caller}-catalog", str(catalogs[caller])]
    args += ["--atarva-catalog-index", str(atarva) + ".tbi"]
    assert invoke(monkeypatch, args) == 0
    report = json.loads((out / "preflight_report.json").read_text())
    assert report["samples"][0]["eligible_callers"] == ["trgt", "longtr", "atarva", "strdust"]
    assert report["samples"][1]["eligible_callers"] == ["longtr", "atarva", "strdust"]


def test_bad_caller_name_is_rejected() -> None:
    try:
        validate_inputs.normalize_callers("trgt,not-a-caller")
    except validate_inputs.ValidationError as error:
        assert "Unknown" in str(error)
    else:
        raise AssertionError("invalid caller accepted")


def test_unsafe_extra_argument_is_rejected(tmp_path: Path) -> None:
    config = caller_config("strdust")
    config["additional_args"] = [";rm"]
    path = tmp_path / "config.yml"; path.write_text(yaml.safe_dump(config))
    try:
        validate_inputs.validate_caller_config("strdust", path, path)
    except validate_inputs.ValidationError as error:
        assert error
    else:
        raise AssertionError("unsafe argument accepted")


def test_caller_config_reference_build_is_optional(tmp_path: Path) -> None:
    config = caller_config("strdust")
    del config["reference_build"]
    path = tmp_path / "config.yml"
    path.write_text(yaml.safe_dump(config))
    validate_inputs.validate_caller_config("strdust", path, path)
