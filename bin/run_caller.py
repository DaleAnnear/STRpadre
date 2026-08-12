#!/usr/bin/env python3
"""Execute one documented tandem-repeat caller with an argv list, never a shell fragment."""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml

PLATFORMS = {"trgt": {"hifi"}, "longtr": {"hifi", "ont"}, "atarva": {"hifi", "ont"}, "strdust": {"hifi", "ont"}}


class CallerError(RuntimeError):
    pass


def load_config(path: Path, caller: str) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text())
    if not isinstance(config, dict) or config.get("caller") != caller:
        raise CallerError(f"Invalid {caller} configuration: {path}")
    if config.get("platforms") != sorted(PLATFORMS[caller]) and set(config.get("platforms", [])) != PLATFORMS[caller]:
        raise CallerError(f"{caller} platform support in configuration has been altered")
    return config


def executable(name: str) -> str:
    value = shutil.which(name)
    if not value:
        raise CallerError(f"Required executable is absent from caller container: {name}")
    return value


def command_output(args: list[str]) -> str:
    completed = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    return completed.stdout.strip()


def run(args: list[str], log: Path, stdout: Path | None = None) -> None:
    with log.open("a") as stderr_handle:
        if stdout:
            with stdout.open("w") as stdout_handle:
                subprocess.run(args, stdout=stdout_handle, stderr=stderr_handle, text=True, check=True)
        else:
            subprocess.run(args, stdout=stderr_handle, stderr=stderr_handle, text=True, check=True)


def karyotype(sex: str, default: str) -> str:
    normalized = sex.lower()
    if normalized in {"male", "m", "xy"}:
        return "XY"
    if normalized in {"female", "f", "xx"}:
        return "XX"
    return default


def option(config: dict[str, Any], name: str, default: Any = None) -> Any:
    return config.get("options", {}).get(name, default)


def append_if(command: list[str], flag: str, value: Any, enabled: bool = True) -> None:
    if enabled and value is not None:
        command.extend([flag, str(value)])


def build_command(caller: str, config: dict[str, Any], ns: argparse.Namespace, outdir: Path) -> tuple[list[str], Path, bool]:
    extra = list(config.get("additional_args", []))
    if caller == "trgt":
        prefix = outdir / "trgt"
        command = ["trgt", "genotype", "--genome", str(ns.reference), "--reads", str(ns.alignment), "--repeats", str(ns.catalog), "--output-prefix", str(prefix), "--sample-name", ns.sample_id, "--karyotype", karyotype(ns.sex, option(config, "karyotype_default", "XX")), "--threads", str(config["threads"])]
        append_if(command, "--preset", option(config, "preset", "wgs"))
        append_if(command, "--genotyper", option(config, "genotyper", "size"))
        append_if(command, "--flank-len", option(config, "flank_len"))
        append_if(command, "--output-flank-len", option(config, "output_flank_len"))
        append_if(command, "--max-depth", option(config, "max_depth"))
        if option(config, "disable_bam_output", False):
            command.append("--disable-bam-output")
        command.extend(extra)
        return command, prefix.with_suffix(".vcf.gz"), False
    if caller == "longtr":
        raw_vcf = outdir / "calls.original.vcf.gz"
        command = ["LongTR", "--bams", str(ns.alignment), "--fasta", str(ns.reference), "--regions", str(ns.catalog), "--tr-vcf", str(raw_vcf)]
        append_if(command, "--min-mapq", option(config, "min_mapq"))
        append_if(command, "--min-mean-qual", option(config, "min_mean_qual"))
        append_if(command, "--max-tr-len", option(config, "max_tr_len"))
        append_if(command, "--min-reads", option(config, "min_reads"))
        append_if(command, "--indel-flank-len", option(config, "indel_flank_len"))
        if config.get("phased_reads_expected") or option(config, "phased_bam", False):
            command.append("--phased-bam")
        if option(config, "output_filters", False):
            command.append("--output-filters")
        haploid = option(config, "haploid_chromosomes", [])
        if haploid:
            append_if(command, "--haploid-chrs", ",".join(haploid))
        if not option(config, "use_lb_tags", False):
            command.append("--lib-from-samp")
        command.extend(extra)
        return command, raw_vcf, False
    if caller == "atarva":
        raw_vcf = outdir / "calls.original.vcf"
        input_format = "cram" if str(ns.alignment).lower().endswith(".cram") else "bam"
        command = ["atarva", "genotype", "--fasta", str(ns.reference), "--bam", str(ns.alignment), "--regions", str(ns.catalog), "--format", input_format, "--vcf", str(raw_vcf), "--threads", str(config["threads"]), "--karyotype", karyotype(ns.sex, option(config, "karyotype_default", "XX"))]
        option_flags = {"map_qual": "--map-qual", "min_reads": "--min-reads", "max_reads": "--max-reads", "snp_dist": "--snp-dist", "snp_count": "--snp-count", "snp_qual": "--snp-qual", "flank": "--flank"}
        for name, flag in option_flags.items():
            append_if(command, flag, option(config, name))
        if config.get("phased_reads_expected"):
            append_if(command, "--haplotag", option(config, "haplotag", "HP"))
        for name, flag in (("decompose", "--decompose"), ("loci_wise", "--loci-wise"), ("amplicon", "--amplicon"), ("somatic", "--somatic")):
            if option(config, name, False):
                command.append(flag)
        command.extend(extra)
        return command, raw_vcf, False
    if caller == "strdust":
        raw_vcf = outdir / "calls.original.vcf"
        command = ["STRdust", "--region-file", str(ns.catalog), "--sample", ns.sample_id, "--threads", str(config["threads"])]
        option_flags = {"minlen": "--minlen", "support": "--support", "consensus_reads": "--consensus-reads", "max_number_reads": "--max-number-reads", "max_locus": "--max-locus"}
        for name, flag in option_flags.items():
            append_if(command, flag, option(config, name))
        if not config.get("phased_reads_expected"):
            command.extend(["--unphased", "--phasing", str(option(config, "phasing", "ward"))])
        if option(config, "find_outliers", False):
            command.append("--find-outliers")
        haploid = option(config, "haploid_chromosomes", [])
        if haploid:
            append_if(command, "--haploid", ",".join(haploid))
        command.extend(extra)
        command.extend([str(ns.reference), str(ns.alignment)])
        return command, raw_vcf, True
    raise CallerError(f"Unsupported caller: {caller}")


def version_probe(caller: str) -> tuple[str, str]:
    probes = {"trgt": ["trgt", "--version"], "longtr": ["LongTR", "--help"], "atarva": ["atarva", "--version"], "strdust": ["STRdust", "--version"]}
    output = command_output(probes[caller])
    return " ".join(probes[caller]), output


LONGTR_DFLANKINDEL_HEADER = (
    '##FORMAT=<ID=DFLANKINDEL,Number=.,Type=String,'
    'Description="LongTR compatibility declaration">\n'
)


def repair_longtr_vcf_header(raw_vcf: Path) -> Path:
    """Add LongTR's omitted DFLANKINDEL header declaration when needed."""
    opener = gzip.open if raw_vcf.suffix == ".gz" else open
    with opener(raw_vcf, "rt") as handle:
        for line in handle:
            if line.startswith("##FORMAT=<ID=DFLANKINDEL,"):
                return raw_vcf
            if line.startswith("#CHROM"):
                break

    repaired_vcf = raw_vcf.with_name(f"{raw_vcf.stem}.with-header{raw_vcf.suffix}")
    with opener(raw_vcf, "rt") as source, opener(repaired_vcf, "wt") as destination:
        for line in source:
            if line.startswith("#CHROM"):
                destination.write(LONGTR_DFLANKINDEL_HEADER)
            destination.write(line)
    return repaired_vcf


def normalise_vcf(raw_vcf: Path, output: Path, log: Path, caller: str) -> None:
    if not raw_vcf.is_file() or raw_vcf.stat().st_size == 0:
        raise CallerError(f"Caller did not create a non-empty VCF: {raw_vcf}")
    source_vcf = repair_longtr_vcf_header(raw_vcf) if caller == "longtr" else raw_vcf
    if source_vcf != raw_vcf:
        with log.open("a") as stderr_handle:
            stderr_handle.write("Added missing LongTR DFLANKINDEL FORMAT header declaration.\n")
    run([executable("bcftools"), "sort", "-Oz", "-o", str(output), str(source_vcf)], log)
    run([executable("tabix"), "-f", "-p", "vcf", str(output)], log)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--caller", choices=PLATFORMS, required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--platform", choices=["hifi", "ont"], required=True)
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--alignment-index", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--sex", default="")
    parser.add_argument("--ploidy", default="")
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    ns = parse_args()
    try:
        if ns.platform not in PLATFORMS[ns.caller]:
            raise CallerError(f"{ns.caller} is not applicable to {ns.platform}")
        for source in (ns.alignment, ns.alignment_index, ns.reference, ns.catalog, ns.config):
            if not source.is_file():
                raise CallerError(f"Required staged input missing: {source}")
        config = load_config(ns.config, ns.caller)
        ns.output_dir.mkdir(parents=True, exist_ok=False)
        command, raw_vcf, stdout_is_vcf = build_command(ns.caller, config, ns, ns.output_dir)
        executable(command[0])
        log = ns.output_dir / "caller.stderr.log"
        probe_command, actual_version = version_probe(ns.caller)
        run(command, log, raw_vcf if stdout_is_vcf else None)
        normalise_vcf(raw_vcf, ns.output_dir / "calls.vcf.gz", log, ns.caller)
        provenance = {"caller": ns.caller, "sample_id": ns.sample_id, "platform": ns.platform, "configured_tool_version": config["tool_version"], "executable_version_probe": probe_command, "executable_version_output": actual_version, "argv": command, "container": config["container"], "native_vcf": raw_vcf.name, "normalized_native_vcf": "calls.vcf.gz"}
        (ns.output_dir / "command.json").write_text(json.dumps(provenance, indent=2) + "\n")
        (ns.output_dir / "versions.yml").write_text(yaml.safe_dump({"caller": ns.caller, "configured_version": config["tool_version"], "executable_probe": probe_command, "executable_output": actual_version, "bcftools": command_output(["bcftools", "--version"]).splitlines()[0], "tabix": command_output(["tabix", "--version"]).splitlines()[0]}, sort_keys=False))
    except (CallerError, OSError, subprocess.CalledProcessError) as exc:
        print(f"CALLER FAILURE [{ns.caller}/{ns.sample_id}]: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
