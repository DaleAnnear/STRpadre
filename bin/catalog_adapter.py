#!/usr/bin/env python3
"""Create documented caller-native repeat catalogs from one canonical manifest.

Only conversions that retain target-caller semantics are made. A canonical primary
motif cannot describe a complex TRGT ``MOTIFS``/``STRUC`` locus, so TRGT output is
limited to simple single-motif loci.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any

try:
    import pysam
except ImportError:
    pysam = None


CALLERS = ("trgt", "longtr", "atarva", "strdust")
REQUIRED_COLUMNS = {
    "locus_id",
    "chr",
    "start",
    "end",
    "coordinate_system",
    "reference_build",
    "primary_motif",
}
DNA = re.compile(r"[ACGTN]+")
NATIVE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.:-]*")
MAPPING_FIELDS = (
    "canonical_locus_id",
    "caller",
    "native_locus_id",
    "chr",
    "start",
    "end",
    "motif",
    "catalog_file",
    "catalog_coordinate_system",
)


class CatalogAdapterError(ValueError):
    """An input catalog cannot be converted without making up semantics."""


def fail(message: str) -> None:
    raise CatalogAdapterError(message)


def normalize_callers(value: str) -> list[str]:
    names = [name.strip().lower() for name in value.split(",") if name.strip()]
    if names == ["all"]:
        return list(CALLERS)
    if not names:
        fail("Caller selection is empty")
    unknown = sorted(set(names) - set(CALLERS))
    if unknown:
        fail(f"Unknown caller(s): {', '.join(unknown)}")
    if len(names) != len(set(names)):
        fail("Duplicate caller names are not allowed")
    return names


def read_fai(path: Path) -> dict[str, int]:
    if not path.is_file():
        fail(f"Reference FAI does not exist: {path}")
    contigs: dict[str, int] = {}
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split("\t")
        if len(fields) < 2 or not fields[0]:
            fail(f"Malformed FAI record at {path}:{line_number}")
        try:
            length = int(fields[1])
        except ValueError:
            fail(f"Invalid FAI length at {path}:{line_number}")
        if length < 1 or fields[0] in contigs:
            fail(f"Invalid or duplicate FAI contig at {path}:{line_number}")
        contigs[fields[0]] = length
    if not contigs:
        fail(f"Reference FAI is empty: {path}")
    return contigs


def load_manifest(path: Path, fai: dict[str, int]) -> tuple[list[dict[str, str]], str]:
    if not path.is_file():
        fail(f"Canonical manifest does not exist: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not REQUIRED_COLUMNS.issubset(reader.fieldnames):
            fail(f"Canonical manifest lacks required columns: {sorted(REQUIRED_COLUMNS)}")
        rows: list[dict[str, str]] = []
        locus_ids: set[str] = set()
        builds: set[str] = set()
        for line_number, raw in enumerate(reader, 2):
            row = {key: (value or "").strip() for key, value in raw.items()}
            locus_id = row["locus_id"]
            if not NATIVE_ID.fullmatch(locus_id) or locus_id in locus_ids:
                fail(f"locus_id must be unique and shell-safe at manifest line {line_number}")
            locus_ids.add(locus_id)
            if row["coordinate_system"] != "0-based-half-open":
                fail("Catalog adaptation requires coordinate_system=0-based-half-open")
            try:
                start, end = int(row["start"]), int(row["end"])
            except ValueError:
                fail(f"Coordinates must be integers at manifest line {line_number}")
            contig = row["chr"]
            if contig not in fai:
                fail(f"Manifest contig {contig!r} at line {line_number} is absent from the FAI")
            if start < 0 or end <= start or end > fai[contig]:
                message = f"Manifest interval {contig}:{start}-{end} at line {line_number} "
                fail(message + "is outside the reference")
            motif = row["primary_motif"].upper()
            if not DNA.fullmatch(motif):
                fail(
                    "primary_motif must be unambiguous DNA/IUPAC N "
                    f"at manifest line {line_number}"
                )
            build = row["reference_build"]
            if not build:
                fail(f"reference_build is empty at manifest line {line_number}")
            for caller in CALLERS:
                configured_id = row.get(f"{caller}_id", "")
                if configured_id and not NATIVE_ID.fullmatch(configured_id):
                    message = f"{caller}_id must be a safe native identifier "
                    fail(message + f"at manifest line {line_number}")
            row["start"] = str(start)
            row["end"] = str(end)
            row["primary_motif"] = motif
            builds.add(build)
            rows.append(row)
    if not rows:
        fail("Canonical manifest has no loci")
    if len(builds) != 1:
        fail("All manifest loci must have the same non-empty reference_build")
    order = {contig: index for index, contig in enumerate(fai)}
    rows.sort(
        key=lambda row: (
            order[row["chr"]],
            int(row["start"]),
            int(row["end"]),
            row["locus_id"],
        )
    )
    return rows, next(iter(builds))


def native_id(row: dict[str, str], caller: str) -> str:
    """Return caller-visible ID or a coordinate identity where the format has no ID."""
    configured = row.get(f"{caller}_id", "")
    if configured:
        return configured
    if caller in {"trgt", "longtr"}:
        return row["locus_id"]
    return f"{row['chr']}:{int(row['start']) + 1}-{row['end']}"


def validate_target_identity(rows: list[dict[str, str]], caller: str) -> None:
    """Reject records that a target format would make indistinguishable."""
    if caller in {"trgt", "longtr"}:
        identities: list[Any] = [native_id(row, caller) for row in rows]
    else:
        identities = [(row["chr"], row["start"], row["end"]) for row in rows]
    if len(identities) != len(set(identities)):
        fail(f"{caller} output would contain duplicate, indistinguishable locus records")


def target_name(caller: str) -> str:
    return {
        "trgt": "trgt.bed",
        "longtr": "longtr.bed",
        "atarva": "atarva.bed.gz",
        "strdust": "strdust.bed",
    }[caller]


def target_coordinate_system(caller: str) -> str:
    return "1-based-closed" if caller == "longtr" else "0-based-half-open"


def write_catalog(
    caller: str,
    rows: list[dict[str, str]],
    directory: Path,
) -> list[dict[str, str]]:
    validate_target_identity(rows, caller)
    if caller == "atarva" and pysam is None:
        fail("pysam is required to bgzip and tabix-index the ATaRVa catalog")
    filename = target_name(caller)
    final_path = directory / filename
    text_path = final_path.with_suffix("") if caller == "atarva" else final_path
    mappings: list[dict[str, str]] = []
    with text_path.open("w", newline="") as handle:
        if caller == "atarva":
            handle.write("#CHROM\tSTART\tEND\tMOTIF\tMOTIF_LEN\n")
        for row in rows:
            contig = row["chr"]
            start = int(row["start"])
            end = int(row["end"])
            motif = row["primary_motif"]
            visible_id = native_id(row, caller)
            if caller == "trgt":
                if row.get("repeat_structure", ""):
                    fail(
                        f"TRGT catalog cannot be generated for complex locus {row['locus_id']}; "
                        "provide a native TRGT catalog with complete MOTIFS and STRUC fields"
                    )
                handle.write(
                    f"{contig}\t{start}\t{end}\tID={visible_id};MOTIFS={motif};STRUC=<TR>\n"
                )
            elif caller == "longtr":
                handle.write(f"{contig}\t{start + 1}\t{end}\t{motif}\t{visible_id}\n")
            elif caller == "atarva":
                handle.write(f"{contig}\t{start}\t{end}\t{motif}\t{len(motif)}\n")
            else:
                handle.write(f"{contig}\t{start}\t{end}\n")
            mappings.append(
                {
                    "canonical_locus_id": row["locus_id"],
                    "caller": caller,
                    "native_locus_id": visible_id,
                    "chr": contig,
                    "start": str(start),
                    "end": str(end),
                    "motif": motif,
                    "catalog_file": filename,
                    "catalog_coordinate_system": target_coordinate_system(caller),
                }
            )
    if caller == "atarva":
        pysam.tabix_compress(str(text_path), str(final_path), force=True)
        pysam.tabix_index(str(final_path), preset="bed", force=True)
        text_path.unlink()
    return mappings


def write_mapping(directory: Path, mappings: list[dict[str, str]]) -> None:
    with (directory / "catalog_locus_mapping.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MAPPING_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(mappings)


def generate(
    rows: list[dict[str, str]],
    callers: list[str],
    reference_build: str,
    source_manifest: Path,
    output_dir: Path,
) -> None:
    if output_dir.exists():
        fail(f"Refusing to overwrite existing output directory: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.tmp-", dir=output_dir.parent)
    )
    try:
        mappings: list[dict[str, str]] = []
        for caller in callers:
            mappings.extend(write_catalog(caller, rows, temporary))
        write_mapping(temporary, mappings)
        trgt_policy = (
            "simple single-motif loci only; complex repeat_structure "
            "requires native catalog"
        )
        report: dict[str, Any] = {
            "schema_version": "1.0",
            "source_manifest": str(source_manifest),
            "reference_build": reference_build,
            "canonical_coordinate_system": "0-based-half-open",
            "locus_count": len(rows),
            "callers": callers,
            "outputs": {
                caller: {
                    "catalog": target_name(caller),
                    "coordinate_system": target_coordinate_system(caller),
                    "locus_count": len(rows),
                    "conversion_policy": (
                        trgt_policy
                        if caller == "trgt"
                        else "explicit documented coordinate and field mapping"
                    ),
                }
                for caller in callers
            },
        }
        report_path = temporary / "catalog_adaptation_report.json"
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, output_dir)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--locus-manifest",
        type=Path,
        required=True,
        help="Canonical TSV locus manifest",
    )
    parser.add_argument(
        "--reference-fai",
        type=Path,
        required=True,
        help="Reference FASTA .fai used to validate contigs and bounds",
    )
    parser.add_argument(
        "--callers",
        default="trgt,longtr,atarva,strdust",
        help="Comma-separated callers or 'all'",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="New directory for generated catalogs and provenance",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        callers = normalize_callers(args.callers)
        fai = read_fai(args.reference_fai)
        rows, reference_build = load_manifest(args.locus_manifest, fai)
        generate(rows, callers, reference_build, args.locus_manifest, args.output_dir)
    except (CatalogAdapterError, OSError) as exc:
        print(f"CATALOG ADAPTER FAILURE: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
