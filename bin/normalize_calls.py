#!/usr/bin/env python3
"""Convert documented caller-specific VCF fields into STRpadre's loss-aware common TSV."""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import yaml

PARSER_SCHEMA_VERSION = "1.0.0"
MISSING = "."
NORMALIZED_COLUMNS = [
    "sample_id", "platform", "canonical_locus_id", "chr", "canonical_start", "canonical_end",
    "reference_build", "motif", "caller", "call_status", "allele1_length_bp", "allele2_length_bp",
    "allele1_copy_number", "allele2_copy_number", "allele1_sequence", "allele2_sequence",
    "genotype_quality", "read_support_allele1", "read_support_allele2", "phase_set",
    "haplotype_labels", "filter_reasons", "native_record_id", "source_vcf", "parser_schema_version",
]


class ParserError(ValueError):
    pass


@dataclass(frozen=True)
class VcfRecord:
    chrom: str
    pos: int
    record_id: str
    ref: str
    alt: list[str]
    filt: str
    info: dict[str, str | bool]
    format: dict[str, str]


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def parse_info(value: str) -> dict[str, str | bool]:
    if value in {"", "."}:
        return {}
    result: dict[str, str | bool] = {}
    for entry in value.split(";"):
        if "=" in entry:
            key, field = entry.split("=", 1)
            result[key] = field
        else:
            result[entry] = True
    return result


def iter_vcf(path: Path, sample_id: str) -> Iterable[VcfRecord]:
    samples: list[str] = []
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                if not samples:
                    raise ParserError(f"{path}: no sample column")
                if sample_id in samples:
                    sample_index = samples.index(sample_id)
                elif len(samples) == 1:
                    sample_index = 0
                else:
                    raise ParserError(f"{path}: requested sample {sample_id} absent from multi-sample VCF")
                continue
            if line.startswith("#"):
                continue
            if not samples:
                raise ParserError(f"{path}: VCF header missing #CHROM")
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                raise ParserError(f"{path}:{line_number}: expected sample VCF record")
            keys = fields[8].split(":")
            values = fields[9 + sample_index].split(":")
            sample_format = {key: values[index] if index < len(values) else MISSING for index, key in enumerate(keys)}
            try:
                pos = int(fields[1])
            except ValueError as exc:
                raise ParserError(f"{path}:{line_number}: non-integer POS") from exc
            yield VcfRecord(fields[0], pos, fields[2], fields[3], [] if fields[4] == "." else fields[4].split(","), fields[6], parse_info(fields[7]), sample_format)


def load_manifest(path: Path) -> tuple[dict[str, dict[str, str]], str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"locus_id", "chr", "start", "end", "coordinate_system", "reference_build", "primary_motif"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ParserError("Manifest lacks required canonical columns")
        loci = {row["locus_id"]: {key: (value or "") for key, value in row.items()} for row in reader}
    builds = {row["reference_build"] for row in loci.values()}
    if len(builds) != 1:
        raise ParserError("Manifest must contain exactly one reference build")
    return loci, next(iter(builds))


def info_text(record: VcfRecord, key: str) -> str:
    value = record.info.get(key, MISSING)
    return MISSING if value is True else str(value)


def native_coordinates(caller: str, record: VcfRecord) -> tuple[int, int, str | None, str | None]:
    if caller == "trgt":
        return record.pos - 1, int(info_text(record, "END")), info_text(record, "TRID"), info_text(record, "MOTIFS").split(",")[0]
    if caller == "longtr":
        return int(info_text(record, "START")) - 1, int(info_text(record, "END")), None if record.record_id == MISSING else record.record_id, info_text(record, "MOTIF").split(",")[0]
    if caller == "atarva":
        return int(info_text(record, "START")), int(info_text(record, "END")), info_text(record, "ID"), info_text(record, "MOTIF")
    if caller == "strdust":
        return record.pos - 1, int(info_text(record, "END")), None, None
    raise ParserError(f"Unknown caller {caller}")


def canon_motif(motif: str) -> str:
    motif = motif.upper()
    rotations = [motif[offset:] + motif[:offset] for offset in range(len(motif))]
    complement = str.maketrans("ACGT", "TGCA")
    reverse_complement = motif.translate(complement)[::-1]
    rotations.extend(reverse_complement[offset:] + reverse_complement[:offset] for offset in range(len(motif)))
    return min(rotations)


def map_locus(caller: str, record: VcfRecord, loci: dict[str, dict[str, str]]) -> dict[str, str]:
    start, end, native_id, motif = native_coordinates(caller, record)
    candidates = []
    if native_id and native_id != MISSING:
        candidates = [row for row in loci.values() if row.get(f"{caller}_id") == native_id or row["locus_id"] == native_id]
    if not candidates:
        candidates = [row for row in loci.values() if row["chr"] == record.chrom and int(row["start"]) == start and int(row["end"]) == end]
    if len(candidates) != 1:
        raise ParserError(f"{caller} VCF record {record.chrom}:{record.pos} cannot be mapped uniquely")
    locus = candidates[0]
    if locus["chr"] != record.chrom or int(locus["start"]) != start or int(locus["end"]) != end:
        raise ParserError(f"{caller} native locus ID does not agree with canonical coordinates")
    if motif and motif != MISSING and canon_motif(locus["primary_motif"]) != canon_motif(motif):
        raise ParserError(f"{caller} native locus does not agree with canonical primary motif")
    return locus


def values(field: str) -> list[str]:
    if field in {"", MISSING}:
        return []
    return re.split(r"[,|/]", field)


def number_values(field: str) -> list[str]:
    result: list[str] = []
    for value in values(field):
        if value in {"", MISSING}:
            result.append(MISSING)
        else:
            try:
                float(value)
            except ValueError:
                result.append(MISSING)
            else:
                result.append(value)
    return result


def pair(items: list[str]) -> tuple[str, str]:
    return (items[0] if items else MISSING, items[1] if len(items) > 1 else MISSING)


def genotype_indices(record: VcfRecord) -> list[str]:
    return values(record.format.get("GT", MISSING))


def sequences(record: VcfRecord, caller: str) -> tuple[str, str]:
    sequence_by_index = [record.ref] + record.alt
    output = []
    for allele in genotype_indices(record)[:2]:
        try:
            index = int(allele)
            sequence = sequence_by_index[index]
        except (ValueError, IndexError):
            sequence = MISSING
        # TRGT documents a leading padding base in REF/ALT. Other callers do not share this convention.
        if caller == "trgt" and sequence != MISSING:
            sequence = sequence[1:]
        output.append(sequence)
    return pair(output)


def copy_numbers_from_length(lengths: tuple[str, str], motif: str) -> tuple[str, str]:
    answer = []
    for length in lengths:
        try:
            number = float(length)
        except ValueError:
            answer.append(MISSING)
            continue
        if number % len(motif) == 0:
            answer.append(str(int(number / len(motif))))
        else:
            answer.append(MISSING)
    return pair(answer)


def parse_record(caller: str, record: VcfRecord, locus: dict[str, str], sample_id: str, platform: str, source_vcf: str) -> dict[str, str]:
    gt = record.format.get("GT", MISSING)
    filtered = record.filt not in {MISSING, "PASS"}
    status = "filtered_call" if filtered else ("valid_no_call" if not genotype_indices(record) else "valid_call")
    if any(item == MISSING for item in genotype_indices(record)):
        status = "filtered_call" if filtered else "valid_no_call"
    motif = locus["primary_motif"]
    sequence1, sequence2 = sequences(record, caller)
    quality = MISSING
    phase_set = record.format.get("PS", MISSING)
    haplotypes = "phased" if "|" in gt else "unphased"
    if caller == "trgt":
        lengths = pair(number_values(record.format.get("AL", MISSING)))
        motif_list = values(info_text(record, "MOTIFS"))
        try:
            motif_index = motif_list.index(motif)
        except ValueError:
            motif_index = -1
        mc = values(record.format.get("MC", MISSING))
        copies = pair([entry.split("_")[motif_index] if motif_index >= 0 and len(entry.split("_")) > motif_index else MISSING for entry in mc])
        supports = pair(number_values(record.format.get("SD", MISSING)))
    elif caller == "longtr":
        differences = pair(number_values(record.format.get("GB", MISSING)))
        base = int(locus["end"]) - int(locus["start"])
        lengths = pair([str(int(base + float(value))) if value != MISSING else MISSING for value in differences])
        copies = copy_numbers_from_length(lengths, motif)
        quality = record.format.get("PQ", record.format.get("Q", MISSING))
        supports = pair(number_values(record.format.get("DP", MISSING)))
    elif caller == "atarva":
        lengths = pair(number_values(record.format.get("AL", MISSING)))
        copies = pair(number_values(record.format.get("CN", MISSING)))
        supports = pair(number_values(record.format.get("SD", MISSING)))
    else:
        lengths = pair(number_values(record.format.get("FRB", MISSING)))
        copies = copy_numbers_from_length(lengths, motif)
        supports = pair(number_values(record.format.get("SUP", MISSING)))
        quality = record.format.get("SC", MISSING)
    return {
        "sample_id": sample_id, "platform": platform, "canonical_locus_id": locus["locus_id"], "chr": locus["chr"],
        "canonical_start": locus["start"], "canonical_end": locus["end"], "reference_build": locus["reference_build"], "motif": motif,
        "caller": caller, "call_status": status, "allele1_length_bp": lengths[0], "allele2_length_bp": lengths[1],
        "allele1_copy_number": copies[0], "allele2_copy_number": copies[1], "allele1_sequence": sequence1, "allele2_sequence": sequence2,
        "genotype_quality": quality, "read_support_allele1": supports[0], "read_support_allele2": supports[1], "phase_set": phase_set,
        "haplotype_labels": haplotypes, "filter_reasons": MISSING if not filtered else record.filt,
        "native_record_id": info_text(record, "TRID") if caller == "trgt" else (info_text(record, "ID") if caller == "atarva" else record.record_id),
        "source_vcf": source_vcf, "parser_schema_version": PARSER_SCHEMA_VERSION,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--caller", choices=["trgt", "longtr", "atarva", "strdust"], required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--platform", choices=["hifi", "ont"], required=True)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--locus-manifest", type=Path, required=True)
    parser.add_argument("--source-vcf", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--versions", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        loci, _ = load_manifest(args.locus_manifest)
        parsed = [parse_record(args.caller, record, map_locus(args.caller, record, loci), args.sample_id, args.platform, args.source_vcf) for record in iter_vcf(args.vcf, args.sample_id)]
        if len({row["canonical_locus_id"] for row in parsed}) != len(parsed):
            raise ParserError("More than one native VCF record mapped to a canonical locus")
        with gzip.open(args.output, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=NORMALIZED_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(parsed)
        args.versions.write_text(yaml.safe_dump({"normalizer": {"caller": args.caller, "parser_schema_version": PARSER_SCHEMA_VERSION, "python": sys.version.split()[0]}}, sort_keys=False))
    except (OSError, ParserError, ValueError) as exc:
        print(f"NORMALIZATION FAILURE [{args.caller}/{args.sample_id}]: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__": raise SystemExit(main())
