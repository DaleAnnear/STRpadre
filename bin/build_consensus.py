#!/usr/bin/env python3
"""Deterministic, provenance-rich consensus for normalized tandem-repeat calls."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import yaml

MISSING = "."
PLATFORMS = {"trgt": {"hifi"}, "longtr": {"hifi", "ont"}, "atarva": {"hifi", "ont"}, "strdust": {"hifi", "ont"}}
REQUIRED_NORMALIZED = {"sample_id", "platform", "canonical_locus_id", "chr", "canonical_start", "canonical_end", "reference_build", "motif", "caller", "call_status", "allele1_length_bp", "allele2_length_bp", "allele1_copy_number", "allele2_copy_number", "read_support_allele1", "read_support_allele2", "source_vcf"}


class ConsensusError(ValueError):
    pass


def as_float(value: str | None) -> float | None:
    if value in {None, "", MISSING}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def display_number(value: float | None) -> str:
    if value is None:
        return MISSING
    return str(int(value)) if value.is_integer() else f"{value:.6g}"


def canonical_motif(motif: str, policy: str) -> str:
    motif = motif.upper()
    rotations = [motif[offset:] + motif[:offset] for offset in range(len(motif))]
    if policy == "exact":
        return motif
    if policy == "rotations_and_reverse_complements":
        complement = str.maketrans("ACGT", "TGCA")
        reverse = motif.translate(complement)[::-1]
        rotations.extend(reverse[offset:] + reverse[:offset] for offset in range(len(reverse)))
    return min(rotations)


def allele_distance(left: float, right: float) -> float:
    """Minimum-cost objective used for deterministic diploid assignment."""
    return abs(left - right)


def assign_alleles(observed: list[float], centers: list[float]) -> tuple[list[tuple[int, int]], float]:
    """Return the minimum-cost bijection; ties retain input order deterministically."""
    if len(observed) != len(centers) or len(observed) not in {1, 2}:
        raise ConsensusError("Only equal-ploidy haploid or diploid allele assignments are supported")
    if len(observed) == 1:
        return [(0, 0)], allele_distance(observed[0], centers[0])
    direct = allele_distance(observed[0], centers[0]) + allele_distance(observed[1], centers[1])
    swap = allele_distance(observed[0], centers[1]) + allele_distance(observed[1], centers[0])
    return ([(0, 0), (1, 1)], direct) if direct <= swap else ([(0, 1), (1, 0)], swap)


def compatible(left: float, right: float, config: dict[str, Any]) -> bool:
    return allele_distance(left, right) <= max(float(config["absolute_allele_length_tolerance_bp"]), float(config["relative_allele_length_tolerance"]) * max(abs(left), abs(right), 1.0))


def median(values: Iterable[float]) -> float:
    return float(statistics.median(list(values)))


def load_manifest(path: Path) -> tuple[dict[str, dict[str, str]], str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = {row["locus_id"]: {key: (value or "") for key, value in row.items()} for row in reader}
    if not rows:
        raise ConsensusError("Canonical locus manifest is empty")
    builds = {row["reference_build"] for row in rows.values()}
    if len(builds) != 1:
        raise ConsensusError("Manifest must use one reference build")
    return rows, next(iter(builds))


def vcf_locus_sort_key(locus: dict[str, str]) -> tuple[int, int, str, int, int, str]:
    """Sort canonical human contigs naturally, then other contigs deterministically."""
    contig = locus["chr"]
    normalized = contig[3:] if contig.lower().startswith("chr") else contig
    if normalized.isdigit() and 1 <= int(normalized) <= 22:
        return (0, int(normalized), "", int(locus["start"]), int(locus["end"]), locus["locus_id"])
    if normalized.upper() == "X":
        return (0, 23, "", int(locus["start"]), int(locus["end"]), locus["locus_id"])
    if normalized.upper() == "Y":
        return (0, 24, "", int(locus["start"]), int(locus["end"]), locus["locus_id"])
    return (1, 0, contig, int(locus["start"]), int(locus["end"]), locus["locus_id"])


def load_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in csv.DictReader(handle)]
    if not rows:
        raise ConsensusError("Sample sheet is empty")
    if len({row["sample_id"] for row in rows}) != len(rows):
        raise ConsensusError("Sample sheet has duplicate sample_id")
    return rows


def read_normalized(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not REQUIRED_NORMALIZED.issubset(reader.fieldnames):
            raise ConsensusError(f"Normalized call file has unsupported schema: {path}")
        return [{key: (value or MISSING) for key, value in row.items()} for row in reader]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def row_alleles(row: dict[str, str]) -> list[float]:
    values = [as_float(row.get("allele1_length_bp")), as_float(row.get("allele2_length_bp"))]
    return [value for value in values if value is not None]


def read_support_ok(row: dict[str, str], alleles: list[float], minimum: int) -> bool:
    if minimum == 0:
        return True
    fields = [as_float(row.get("read_support_allele1")), as_float(row.get("read_support_allele2"))]
    return all(value is not None and value >= minimum for value in fields[:len(alleles)])


def eligible_record(row: dict[str, str], config: dict[str, Any]) -> tuple[bool, str]:
    if row["call_status"] == "valid_call":
        status = "candidate"
    elif row["call_status"] == "filtered_call" and config["filtered_calls_may_contribute"]:
        status = "candidate_filtered"
    else:
        return False, row["call_status"]
    alleles = row_alleles(row)
    if not alleles:
        return False, "valid_no_call"
    if not read_support_ok(row, alleles, int(config["minimum_read_support"])):
        return False, "low_read_support"
    return True, status


def consensus_for_rows(rows: list[dict[str, str]], config: dict[str, Any]) -> dict[str, Any]:
    """Harmonize a sample/locus call set with explicit assignment and robust medians."""
    candidates: list[dict[str, str]] = []
    reasons: list[str] = []
    for row in sorted(rows, key=lambda item: (item["caller"], item.get("source_vcf", ""))):
        contributes, reason = eligible_record(row, config)
        row["_consensus_reason"] = reason
        if contributes:
            candidates.append(row)
        elif reason not in {"valid_no_call", "filtered_call"}:
            reasons.append(f"{row['caller']}:{reason}")
    if not candidates:
        return {"status": "no_consensus", "alleles": [], "support": [], "ranges": [], "eligible": 0, "called": 0, "reasons": sorted(set(reasons)), "contributors": []}
    ploidies = {len(row_alleles(row)) for row in candidates}
    if len(ploidies) != 1 and config["haploid_diploid_policy"] == "require_same_ploidy":
        return {"status": "ploidy_mismatch", "alleles": [], "support": [], "ranges": [], "eligible": len(candidates), "called": len(candidates), "reasons": ["haploid_diploid_mismatch"], "contributors": candidates}
    target_ploidy = max(ploidies)
    candidates = [row for row in candidates if len(row_alleles(row)) == target_ploidy]
    clusters: list[list[tuple[str, float, dict[str, str]]]] = [[(candidates[0]["caller"], value, candidates[0])] for value in row_alleles(candidates[0])]
    assigned = {id(candidates[0])}
    for row in candidates[1:]:
        observed = row_alleles(row)
        centers = [median(value for _, value, _ in cluster) for cluster in clusters]
        pairs, _ = assign_alleles(observed, centers)
        all_compatible = True
        for observed_index, center_index in pairs:
            if not compatible(observed[observed_index], centers[center_index], config):
                all_compatible = False
                reasons.append(f"{row['caller']}:allele_outside_tolerance")
        if all_compatible:
            for observed_index, center_index in pairs:
                clusters[center_index].append((row["caller"], observed[observed_index], row))
            assigned.add(id(row))
    centers = [median(value for _, value, _ in cluster) for cluster in clusters]
    supports = [len({caller for caller, _, _ in cluster}) for cluster in clusters]
    ranges = [(min(value for _, value, _ in cluster), max(value for _, value, _ in cluster)) for cluster in clusters]
    status = "concordant" if all(support >= int(config["minimum_caller_support"]) for support in supports) and not reasons else "discordant"
    if any(support < int(config["minimum_caller_support"]) for support in supports):
        status = "insufficient_support" if not reasons else "discordant"
    return {"status": status, "alleles": centers, "support": supports, "ranges": ranges, "eligible": len(candidates), "called": len(assigned), "reasons": sorted(set(reasons)), "contributors": candidates}


def consensus_copy_numbers(rows: list[dict[str, str]], motif: str, allele_count: int) -> list[str]:
    columns = ("allele1_copy_number", "allele2_copy_number")
    result = []
    for index in range(allele_count):
        numbers = [as_float(row.get(columns[index])) for row in rows]
        numbers = [number for number in numbers if number is not None]
        result.append(display_number(median(numbers)) if numbers else MISSING)
    return result


def write_gzip_tsv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def vcf_header(samples: list[str]) -> str:
    return "\n".join([
        "##fileformat=VCFv4.3",
        "##ALT=<ID=STR,Description=\"Tandem repeat; allele sizes are in FORMAT/AL rather than ALT because callers do not share a sequence representation\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Canonical 0-based half-open end coordinate represented as VCF END\">",
        "##INFO=<ID=MOTIF,Number=1,Type=String,Description=\"Canonical primary motif\">",
        "##INFO=<ID=CALLERS,Number=.,Type=String,Description=\"Selected callers contributing compatible calls\">",
        "##INFO=<ID=NELIG,Number=1,Type=Integer,Description=\"Eligible contributing calls\">",
        "##INFO=<ID=NCALL,Number=1,Type=Integer,Description=\"Calls assigned to consensus clusters\">",
        "##INFO=<ID=CSTATUS,Number=1,Type=String,Description=\"Consensus concordance status\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Intentionally no-call: a universal ALT-indexed genotype is not scientifically valid across callers; use AL and CN\">",
        "##FORMAT=<ID=AL,Number=.,Type=Integer,Description=\"Consensus repeat allele lengths in bp\">",
        "##FORMAT=<ID=CN,Number=.,Type=String,Description=\"Consensus copy numbers when callers provide comparable copy-number units\">",
        "##FORMAT=<ID=SC,Number=.,Type=Integer,Description=\"Caller support count for each consensus allele\">",
        "##FORMAT=<ID=DS,Number=.,Type=String,Description=\"Per-allele min-max caller length range in bp\">",
        "##FORMAT=<ID=CS,Number=1,Type=String,Description=\"Consensus status\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" + ("\t" + "\t".join(samples) if samples else ""),
    ]) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--normalized-files", nargs="*", type=Path, default=[])
    parser.add_argument("--locus-manifest", type=Path, required=True)
    parser.add_argument("--samplesheet", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--callers", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        config = yaml.safe_load(args.config.read_text())
        callers = [value.strip() for value in args.callers.split(",") if value.strip()]
        if not callers or len(callers) != len(set(callers)) or set(callers) - set(PLATFORMS):
            raise ConsensusError("Invalid selected caller list")
        loci, build = load_manifest(args.locus_manifest)
        samples = load_samples(args.samplesheet)
        rows = [row for path in args.normalized_files for row in read_normalized(path)]
        seen = set()
        for row in rows:
            key = (row["sample_id"], row["caller"], row["canonical_locus_id"])
            if key in seen:
                raise ConsensusError(f"Duplicate normalized call for {key}")
            seen.add(key)
            if row["canonical_locus_id"] not in loci or row["reference_build"] != build:
                raise ConsensusError(f"Call lacks canonical build/locus compatibility: {key}")
            locus = loci[row["canonical_locus_id"]]
            if row["chr"] != locus["chr"] or int(row["canonical_start"]) != int(locus["start"]) or int(row["canonical_end"]) != int(locus["end"]):
                raise ConsensusError(f"Coordinate mismatch for {key}")
            if canonical_motif(row["motif"], config["motif_equivalence"]) != canonical_motif(locus["primary_motif"], config["motif_equivalence"]):
                raise ConsensusError(f"Motif mismatch for {key}")
        by_key = {(row["sample_id"], row["caller"], row["canonical_locus_id"]): row for row in rows}
        matrix_rows, consensus_rows, discordant_rows = [], [], []
        summaries: dict[str, Counter] = defaultdict(Counter)
        per_sample_locus: dict[tuple[str, str], dict[str, Any]] = {}
        for sample in samples:
            for locus_id, locus in loci.items():
                selected_rows = []
                caller_statuses = []
                for caller in callers:
                    key = (sample["sample_id"], caller, locus_id)
                    row = by_key.get(key)
                    if caller not in callers:
                        status, contributes = "not_selected", False
                    elif sample["platform"] not in PLATFORMS[caller]:
                        status, contributes = "not_applicable_to_platform", False
                    elif row is None:
                        status, contributes = "missing_locus", False
                    else:
                        contributes, status = eligible_record(row, config)
                        selected_rows.append(row)
                    matrix_rows.append({"sample_id": sample["sample_id"], "platform": sample["platform"], "canonical_locus_id": locus_id, "caller": caller, "caller_status": status, "contributes": str(contributes).lower(), "allele1_length_bp": row["allele1_length_bp"] if row else MISSING, "allele2_length_bp": row["allele2_length_bp"] if row else MISSING, "source_vcf": row["source_vcf"] if row else MISSING})
                    caller_statuses.append(status)
                outcome = consensus_for_rows(selected_rows, config)
                copies = consensus_copy_numbers(outcome["contributors"], locus["primary_motif"], len(outcome["alleles"]))
                ranges = [f"{display_number(low)}-{display_number(high)}" for low, high in outcome["ranges"]]
                row = {"sample_id": sample["sample_id"], "platform": sample["platform"], "canonical_locus_id": locus_id, "chr": locus["chr"], "start": locus["start"], "end": locus["end"], "motif": locus["primary_motif"], "consensus_status": outcome["status"], "consensus_allele1_length_bp": display_number(outcome["alleles"][0]) if outcome["alleles"] else MISSING, "consensus_allele2_length_bp": display_number(outcome["alleles"][1]) if len(outcome["alleles"]) > 1 else MISSING, "consensus_allele1_copy_number": copies[0] if copies else MISSING, "consensus_allele2_copy_number": copies[1] if len(copies) > 1 else MISSING, "allele1_supporting_callers": str(outcome["support"][0]) if outcome["support"] else MISSING, "allele2_supporting_callers": str(outcome["support"][1]) if len(outcome["support"]) > 1 else MISSING, "allele1_length_range_bp": ranges[0] if ranges else MISSING, "allele2_length_range_bp": ranges[1] if len(ranges) > 1 else MISSING, "number_eligible_callers": str(outcome["eligible"]), "number_callers_producing_call": str(outcome["called"]), "contributing_callers": ",".join(sorted({item["caller"] for item in outcome["contributors"]})) or MISSING, "discordance_reasons": ";".join(outcome["reasons"]) or MISSING, "caller_statuses": ";".join(f"{caller}:{status}" for caller, status in zip(callers, caller_statuses))}
                consensus_rows.append(row)
                per_sample_locus[(sample["sample_id"], locus_id)] = row
                summaries[locus_id][outcome["status"]] += 1
                if outcome["status"] not in {"concordant", "no_consensus"}:
                    discordant_rows.append(row)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        consensus_columns = list(consensus_rows[0]) if consensus_rows else ["sample_id", "canonical_locus_id", "consensus_status"]
        write_gzip_tsv(args.output_dir / "consensus.tsv.gz", consensus_columns, consensus_rows)
        write_gzip_tsv(args.output_dir / "caller_matrix.tsv.gz", ["sample_id", "platform", "canonical_locus_id", "caller", "caller_status", "contributes", "allele1_length_bp", "allele2_length_bp", "source_vcf"], matrix_rows)
        write_gzip_tsv(args.output_dir / "discordant_calls.tsv.gz", consensus_columns, discordant_rows)
        locus_summary_rows = [{"canonical_locus_id": locus_id, "chr": loci[locus_id]["chr"], "start": loci[locus_id]["start"], "end": loci[locus_id]["end"], "motif": loci[locus_id]["primary_motif"], "concordant": count["concordant"], "discordant": count["discordant"], "insufficient_support": count["insufficient_support"], "no_consensus": count["no_consensus"], "ploidy_mismatch": count["ploidy_mismatch"]} for locus_id, count in sorted(summaries.items())]
        write_gzip_tsv(args.output_dir / "locus_summary.tsv.gz", list(locus_summary_rows[0]) if locus_summary_rows else ["canonical_locus_id"], locus_summary_rows)
        with (args.output_dir / "consensus.vcf").open("w") as vcf:
            sample_ids = [sample["sample_id"] for sample in samples]
            vcf.write(vcf_header(sample_ids))
            for locus_id, locus in sorted(loci.items(), key=lambda item: vcf_locus_sort_key(item[1])):
                sample_values, statuses, contributors, eligible, called = [], [], set(), 0, 0
                for sample in samples:
                    result = per_sample_locus[(sample["sample_id"], locus_id)]
                    alleles = [result["consensus_allele1_length_bp"]] + ([] if result["consensus_allele2_length_bp"] == MISSING else [result["consensus_allele2_length_bp"]])
                    copies = [result["consensus_allele1_copy_number"]] + ([] if result["consensus_allele2_copy_number"] == MISSING else [result["consensus_allele2_copy_number"]])
                    supports = [result["allele1_supporting_callers"]] + ([] if result["allele2_supporting_callers"] == MISSING else [result["allele2_supporting_callers"]])
                    ranges = [result["allele1_length_range_bp"]] + ([] if result["allele2_length_range_bp"] == MISSING else [result["allele2_length_range_bp"]])
                    sample_values.append(":".join(["./.", ",".join(alleles), ",".join(copies), ",".join(supports), ",".join(ranges), result["consensus_status"]]))
                    statuses.append(result["consensus_status"])
                    contributors.update([] if result["contributing_callers"] == MISSING else result["contributing_callers"].split(","))
                    eligible += int(result["number_eligible_callers"])
                    called += int(result["number_callers_producing_call"])
                summary_status = "concordant" if set(statuses) == {"concordant"} else "mixed"
                info = f"END={locus['end']};MOTIF={locus['primary_motif']};CALLERS={','.join(sorted(contributors)) or '.'};NELIG={eligible};NCALL={called};CSTATUS={summary_status}"
                vcf.write("\t".join([locus["chr"], str(int(locus["start"]) + 1), locus_id, "N", "<STR>", ".", "PASS", info, "GT:AL:CN:SC:DS:CS", *sample_values]) + "\n")
        provenance = {"reference_build": build, "selected_callers": callers, "consensus_config": config, "normalization_inputs": [{"path": str(path), "sha256": sha256(path)} for path in args.normalized_files], "algorithm": {"diploid_assignment": "minimum total absolute allele-length difference with direct-order tie break", "statistic": config["consensus_statistic"], "sequence_consensus": "disabled; caller sequence representations are not declared comparable"}}
        (args.output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
        status_counts = Counter(row["consensus_status"] for row in consensus_rows)
        (args.output_dir / "qc_report.md").write_text("# STRpadre consensus QC\n\n" + f"- Samples: {len(samples)}\n- Canonical loci: {len(loci)}\n- Selected callers: {', '.join(callers)}\n- Consensus statuses: {dict(sorted(status_counts.items()))}\n- Sequence-level consensus: disabled; see companion TSV for caller-specific sequences.\n")
        (args.output_dir / "versions.yml").write_text(yaml.safe_dump({"consensus": {"script": "build_consensus.py", "algorithm_version": "1.0.0", "python": sys.version.split()[0]}}, sort_keys=False))
    except (OSError, KeyError, TypeError, ValueError, ConsensusError) as exc:
        print(f"CONSENSUS FAILURE: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
