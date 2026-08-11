#!/usr/bin/env python3
"""Fail-fast, schema-driven preflight for the STRpadre workflow."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import sys
from pathlib import Path
from typing import Any, Iterable

import jsonschema
import yaml

try:
    import pysam
except ImportError:  # Unit tests may exercise the schema layer without pysam.
    pysam = None

CALLERS = ("trgt", "longtr", "atarva", "strdust")
PLATFORMS = {"trgt": {"hifi"}, "longtr": {"hifi", "ont"}, "atarva": {"hifi", "ont"}, "strdust": {"hifi", "ont"}}
OPTION_KEYS = {
    "trgt": {"preset", "genotyper", "flank_len", "output_flank_len", "max_depth", "disable_bam_output", "karyotype_default"},
    "longtr": {"min_mapq", "min_mean_qual", "max_tr_len", "min_reads", "indel_flank_len", "phased_bam", "output_filters", "haploid_chromosomes", "use_lb_tags"},
    "atarva": {"map_qual", "min_reads", "max_reads", "snp_dist", "snp_count", "snp_qual", "flank", "haplotag", "decompose", "loci_wise", "amplicon", "somatic", "karyotype_default"},
    "strdust": {"minlen", "support", "consensus_reads", "max_number_reads", "max_locus", "find_outliers", "phasing", "haploid_chromosomes"},
}
SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


class ValidationError(ValueError):
    pass


def fail(message: str) -> None:
    raise ValidationError(message)


def load_yaml(path: Path) -> dict[str, Any]:
    try:
        data = yaml.safe_load(path.read_text())
    except (OSError, yaml.YAMLError) as exc:
        fail(f"Cannot read YAML {path}: {exc}")
    if not isinstance(data, dict):
        fail(f"YAML document must be an object: {path}")
    return data


def validate_schema(data: dict[str, Any], schema_name: str, source: Path) -> None:
    schema_path = Path(__file__).resolve().parents[1] / "schemas" / schema_name
    try:
        jsonschema.Draft202012Validator(json.loads(schema_path.read_text())).validate(data)
    except (OSError, json.JSONDecodeError, jsonschema.ValidationError) as exc:
        fail(f"Invalid {source}: {exc.message if hasattr(exc, 'message') else exc}")


def normalize_callers(value: str) -> list[str]:
    callers = [item.strip().lower() for item in value.split(",") if item.strip()]
    if not callers:
        fail("Caller selection is empty")
    unknown = sorted(set(callers) - set(CALLERS))
    if unknown:
        fail(f"Unknown caller(s): {', '.join(unknown)}")
    if len(callers) != len(set(callers)):
        fail("Duplicate caller names are not allowed")
    return callers


def read_fai(path: Path) -> dict[str, int]:
    contigs: dict[str, int] = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split("\t")
        if len(fields) < 2:
            fail(f"Malformed FAI {path}:{number}")
        if fields[0] in contigs:
            fail(f"Duplicate FAI contig {fields[0]}")
        try:
            contigs[fields[0]] = int(fields[1])
        except ValueError:
            fail(f"Invalid FAI length at {path}:{number}")
    if not contigs:
        fail(f"Reference FAI is empty: {path}")
    return contigs


def read_manifest(path: Path) -> tuple[dict[str, dict[str, str]], str]:
    required = {"locus_id", "chr", "start", "end", "coordinate_system", "reference_build", "primary_motif"}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            fail(f"Canonical locus manifest lacks required columns: {sorted(required)}")
        loci: dict[str, dict[str, str]] = {}
        builds: set[str] = set()
        for line, row in enumerate(reader, 2):
            locus_id = (row.get("locus_id") or "").strip()
            if not locus_id or locus_id in loci:
                fail(f"Manifest locus_id is empty or duplicated at line {line}")
            try:
                start, end = int(row["start"]), int(row["end"])
            except ValueError:
                fail(f"Manifest coordinates must be integers at line {line}")
            if start < 0 or end <= start:
                fail(f"Manifest has invalid half-open interval at line {line}")
            if row["coordinate_system"] != "0-based-half-open":
                fail("Initial implementation requires canonical manifest coordinate_system=0-based-half-open")
            motif = (row["primary_motif"] or "").upper()
            if not re.fullmatch(r"[ACGTN]+", motif):
                fail(f"Manifest primary_motif must be DNA at line {line}")
            row = {key: (value or "").strip() for key, value in row.items()}
            row["start"], row["end"], row["primary_motif"] = str(start), str(end), motif
            loci[locus_id] = row
            builds.add(row["reference_build"])
    if len(builds) != 1 or not next(iter(builds), ""):
        fail("All manifest loci must share one non-empty reference_build")
    return loci, next(iter(builds))


def valid_additional_args(caller: str, args: list[str]) -> None:
    forbidden = {
        "trgt": {"-g", "--genome", "-r", "--reads", "-b", "--repeats", "-o", "--output-prefix"},
        "longtr": {"--bams", "--fasta", "--regions", "--tr-vcf"},
        "atarva": {"-f", "--fasta", "-b", "--bam", "-r", "--regions", "-o", "--vcf"},
        "strdust": {"-r", "--region", "-R", "--region-file", "--pathogenic", "--sample"},
    }[caller]
    for token in args:
        if not isinstance(token, str) or not re.fullmatch(r"[^\s;&|`$<>\\]+", token):
            fail(f"{caller}: additional_args must be individual safe argv tokens")
        if token.split("=", 1)[0] in forbidden:
            fail(f"{caller}: additional_args may not replace workflow-owned input/output argument {token}")


def validate_caller_config(caller: str, config_path: Path, expected_catalog: Path) -> dict[str, Any]:
    cfg = load_yaml(config_path)
    validate_schema(cfg, "caller.schema.json", config_path)
    if cfg["caller"] != caller:
        fail(f"{config_path}: caller must be {caller}")
    if set(cfg["platforms"]) != PLATFORMS[caller]:
        fail(f"{caller}: platforms must be exactly {sorted(PLATFORMS[caller])}; platform support is not user-overridable")
    if cfg["threads"] != 1 and caller == "longtr":
        fail("LongTR does not support multi-threaded calling; configure threads: 1 and parallelize samples/loci externally")
    if ":latest" in cfg["container"] or "@sha256" not in cfg["container"] and ":" not in cfg["container"]:
        fail(f"{caller}: container must carry a non-latest immutable version tag or digest")
    if cfg["tool_version"] not in cfg["container"]:
        fail(f"{caller}: container tag must include configured tool_version {cfg['tool_version']}")
    if not expected_catalog.is_file():
        fail(f"{caller}: staged native catalog is missing: {expected_catalog}")
    unknown_options = set(cfg["options"]) - OPTION_KEYS[caller]
    if unknown_options:
        fail(f"{caller}: unknown options are forbidden: {sorted(unknown_options)}")
    if caller == "longtr" and not isinstance(cfg["options"].get("use_lb_tags", False), bool):
        fail("longtr: options.use_lb_tags must be true or false")
    valid_additional_args(caller, cfg["additional_args"])
    return cfg


def canonical_motif(motif: str) -> str:
    motif = motif.upper()
    return min(motif[index:] + motif[:index] for index in range(len(motif)))


def find_locus(loci: dict[str, dict[str, str]], caller: str, native_id: str | None, contig: str, start: int, end: int, motif: str | None) -> dict[str, str]:
    if native_id:
        matches = [row for row in loci.values() if row.get(f"{caller}_id") == native_id or row["locus_id"] == native_id]
        if len(matches) == 1:
            row = matches[0]
            if row["chr"] != contig or int(row["start"]) != start or int(row["end"]) != end:
                fail(f"{caller} native ID {native_id} has incompatible coordinates")
            if motif and canonical_motif(row["primary_motif"]) != canonical_motif(motif):
                fail(f"{caller} native ID {native_id} has incompatible primary motif")
            return row
        if len(matches) > 1:
            fail(f"{caller} native ID {native_id} maps to multiple canonical loci")
    matches = [row for row in loci.values() if row["chr"] == contig and int(row["start"]) == start and int(row["end"]) == end]
    if motif:
        matches = [row for row in matches if canonical_motif(row["primary_motif"]) == canonical_motif(motif)]
    if len(matches) != 1:
        fail(f"{caller} catalog locus {contig}:{start}-{end} cannot be mapped uniquely to manifest")
    return matches[0]


def noncomment_lines(path: Path) -> Iterable[list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                yield line.rstrip("\n").split("\t")


def catalog_records(caller: str, path: Path, loci: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for fields in noncomment_lines(path):
        if caller == "trgt":
            if len(fields) < 4:
                fail("TRGT catalog requires BED columns plus structured fourth column")
            metadata = dict(item.split("=", 1) for item in fields[3].split(";") if "=" in item)
            if not {"ID", "MOTIFS", "STRUC"}.issubset(metadata):
                fail("TRGT catalog fourth column requires ID, MOTIFS, and STRUC")
            contig, start, end, native_id, motif = fields[0], int(fields[1]), int(fields[2]), metadata["ID"], metadata["MOTIFS"].split(",")[0]
        elif caller == "longtr":
            if len(fields) < 4:
                fail("LongTR catalog requires at least four columns")
            contig, start, end, motif = fields[0], int(fields[1]) - 1, int(fields[2]), fields[3].split(",")[0]
            native_id = fields[4] if len(fields) > 4 else None
        elif caller == "atarva":
            if len(fields) < 5:
                fail("ATaRVa catalog requires chromosome, start, end, motif, and motif length")
            contig, start, end, motif, native_id = fields[0], int(fields[1]), int(fields[2]), fields[3], None
            if int(fields[4]) != len(motif):
                fail(f"ATaRVa catalog motif length mismatch at {contig}:{start}-{end}")
        else:
            if len(fields) < 3:
                fail("STRdust catalog requires BED columns chrom, start, end")
            contig, start, end, native_id, motif = fields[0], int(fields[1]), int(fields[2]), None, None
        row = find_locus(loci, caller, native_id, contig, start, end, motif)
        records.append({"canonical_locus_id": row["locus_id"], "caller": caller, "native_locus_id": native_id or f"{contig}:{start + 1}-{end}", "chr": contig, "start": str(start), "end": str(end), "motif": motif or row["primary_motif"]})
    if not records:
        fail(f"{caller} native catalog contains no records")
    if len({record["canonical_locus_id"] for record in records}) != len(records):
        fail(f"{caller} native catalog maps multiple rows to one canonical locus")
    return records


def stage_path(original: str, staged_dir: Path | None) -> Path:
    path = Path(original)
    if staged_dir and not path.is_file():
        candidate = staged_dir / path.name
        if candidate.is_file():
            return candidate
    return path


def validate_samples(path: Path, fai: dict[str, int], required_contigs: set[str], reference: Path, staged_dir: Path | None, check_headers: bool) -> list[dict[str, str]]:
    required = {"sample_id", "platform", "alignment", "alignment_index"}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            fail(f"Sample sheet requires columns: {sorted(required)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        fail("Sample sheet is empty")
    names = set()
    basenames = set()
    for row in rows:
        sample = row["sample_id"]
        if not SAFE_SAMPLE_ID.fullmatch(sample) or sample in names:
            fail(f"sample_id must be unique and shell-safe: {sample!r}")
        names.add(sample)
        if row["platform"] not in {"hifi", "ont"}:
            fail(f"{sample}: platform must be hifi or ont")
        if row.get("ploidy") and row["ploidy"] not in {"1", "2"}:
            fail(f"{sample}: ploidy must be 1 or 2 when supplied")
        for key in ("alignment", "alignment_index"):
            item = stage_path(row[key], staged_dir)
            if not item.is_file():
                fail(f"{sample}: {key} does not exist or was not staged: {row[key]}")
            if item.name in basenames:
                fail(f"Alignment/index basenames must be unique for preflight staging: {item.name}")
            basenames.add(item.name)
            row[f"_{key}_resolved"] = str(item)
        if check_headers:
            if pysam is None:
                fail("pysam is required for alignment-header validation")
            try:
                with pysam.AlignmentFile(row["_alignment_resolved"], reference_filename=str(reference)) as alignment:
                    header = alignment.header.to_dict()
            except (OSError, ValueError) as exc:
                fail(f"{sample}: cannot open coordinate-sorted alignment: {exc}")
            sq = {entry.get("SN"): entry.get("LN") for entry in header.get("SQ", [])}
            if not sq:
                fail(f"{sample}: alignment has no @SQ records")
            missing = sorted(required_contigs - set(sq))
            if missing:
                fail(f'{sample}: alignment lacks contigs required by the locus manifest: {", ".join(missing)}')
            mismatches = [contig for contig in required_contigs if sq[contig] != fai[contig]]
            if mismatches:
                fail(f"{sample}: alignment and reference disagree on contig lengths: {mismatches[:3]}")
            if (header.get("HD") or {}).get("SO") != "coordinate":
                fail(f"{sample}: alignment @HD SO must be coordinate")
    return rows


def write_outputs(output_dir: Path, mappings: list[dict[str, str]], report: dict[str, Any]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    columns = ["canonical_locus_id", "caller", "native_locus_id", "chr", "start", "end", "motif"]
    with (output_dir / "locus_mapping.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(mappings, key=lambda record: (record["caller"], record["canonical_locus_id"])))
    (output_dir / "preflight_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (output_dir / "versions.yml").write_text(yaml.safe_dump({"validation": {"script": "validate_inputs.py", "schema_version": "1.0.0", "python": sys.version.split()[0]}}, sort_keys=False))


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("samplesheet", "reference", "reference-fai", "locus-manifest", "trgt-config", "trgt-catalog", "longtr-config", "longtr-catalog", "atarva-config", "atarva-catalog", "atarva-catalog-index", "strdust-config", "strdust-catalog", "consensus-config"):
        parser.add_argument(f"--{key}", required=True)
    parser.add_argument("--callers", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--staged-alignment-dir", type=Path)
    parser.add_argument("--skip-alignment-header-check", action="store_true")
    return parser


def main() -> int:
    args = make_parser().parse_args()
    try:
        callers = normalize_callers(args.callers)
        reference, fai_path = Path(args.reference), Path(args.reference_fai)
        if not reference.is_file() or not fai_path.is_file():
            fail("Reference FASTA and .fai must exist")
        fai = read_fai(fai_path)
        loci, build = read_manifest(Path(args.locus_manifest))
        required_contigs = {row["chr"] for row in loci.values()}
        missing_reference_contigs = sorted(required_contigs - set(fai))
        if missing_reference_contigs:
            fail(f'Reference FAI lacks contigs required by the locus manifest: {", ".join(missing_reference_contigs)}')
        samples = validate_samples(Path(args.samplesheet), fai, required_contigs, reference, args.staged_alignment_dir, not args.skip_alignment_header_check)
        configs = {}
        mappings: list[dict[str, str]] = []
        for caller in CALLERS:
            config = validate_caller_config(caller, Path(getattr(args, f"{caller}_config")), Path(getattr(args, f"{caller}_catalog")))
            configs[caller] = config
            if caller == "atarva":
                catalog_index = Path(args.atarva_catalog_index)
                if not catalog_index.is_file() or catalog_index.suffix != ".tbi":
                    fail("ATaRVa requires a bgzip-compressed catalog with a staged .tbi index")
            mappings.extend(catalog_records(caller, Path(getattr(args, f"{caller}_catalog")), loci))
        consensus = load_yaml(Path(args.consensus_config))
        validate_schema(consensus, "consensus.schema.json", Path(args.consensus_config))
        eligible = {row["sample_id"]: [caller for caller in callers if row["platform"] in PLATFORMS[caller]] for row in samples}
        write_outputs(Path(args.output_dir), mappings, {"status": "valid", "callers": callers, "reference_build": build, "samples": [{"sample_id": row["sample_id"], "platform": row["platform"], "eligible_callers": eligible[row["sample_id"]]} for row in samples], "catalog_locus_counts": {caller: len([mapping for mapping in mappings if mapping["caller"] == caller]) for caller in CALLERS}})
    except ValidationError as exc:
        print(f"INPUT VALIDATION ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
