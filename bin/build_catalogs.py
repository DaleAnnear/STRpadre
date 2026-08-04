#!/usr/bin/env python3
"""Generate only documented native catalog formats from the canonical manifest."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

try:
    import pysam
except ImportError:
    pysam = None


def load(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = [{key: (value or "") for key, value in row.items()} for row in csv.DictReader(handle, delimiter="\t")]
    required = {"locus_id", "chr", "start", "end", "coordinate_system", "primary_motif"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Canonical manifest lacks required fields")
    for row in rows:
        if row["coordinate_system"] != "0-based-half-open":
            raise ValueError("Only 0-based-half-open canonical manifests are supported")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--caller", choices=["longtr", "atarva", "strdust", "trgt"], required=True)
    parser.add_argument("--locus-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-simple-trgt", action="store_true", help="Allow TRGT conversion only for loci without repeat_structure")
    args = parser.parse_args()
    try:
        rows = load(args.locus_manifest)
        plain = args.output.with_suffix("") if args.caller == "atarva" and args.output.suffix == ".gz" else args.output
        with plain.open("w") as handle:
            if args.caller == "atarva":
                handle.write("#CHROM\tSTART\tEND\tMOTIF\tMOTIF_LEN\n")
            for row in rows:
                start, end, motif = int(row["start"]), int(row["end"]), row["primary_motif"]
                native_id = row.get(f"{args.caller}_id") or row["locus_id"]
                if args.caller == "longtr":
                    handle.write(f"{row['chr']}\t{start + 1}\t{end}\t{motif}\t{native_id}\n")
                elif args.caller == "atarva":
                    handle.write(f"{row['chr']}\t{start}\t{end}\t{motif}\t{len(motif)}\n")
                elif args.caller == "strdust":
                    handle.write(f"{row['chr']}\t{start}\t{end}\n")
                else:
                    if not args.allow_simple_trgt or row.get("repeat_structure"):
                        raise ValueError("TRGT conversion is limited to simple loci without repeat_structure; provide a native TRGT catalog otherwise")
                    handle.write(f"{row['chr']}\t{start}\t{end}\tID={native_id};MOTIFS={motif};STRUC=<TR>\n")
        if args.caller == "atarva":
            if pysam is None:
                raise ValueError("pysam is required to bgzip and tabix-index an ATaRVa catalog")
            pysam.tabix_compress(str(plain), str(args.output), force=True)
            pysam.tabix_index(str(args.output), preset="bed", force=True)
            plain.unlink()
    except (OSError, ValueError) as exc:
        print(f"CATALOG ADAPTER FAILURE: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
