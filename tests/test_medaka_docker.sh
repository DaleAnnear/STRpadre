#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
docker run --rm -i -v "$root_dir:/pipeline:ro" -w /tmp local/strpadre-medaka:2.2.2 bash -eu <<'BASH'
python3 - <<'PY'
import random
from pathlib import Path
import pysam

rng = random.Random(43)
left, right = ("".join(rng.choices("ACGT", k=200)) for _ in range(2))
reference = left + "TAC" + "CAG" * 10 + "CCG" + right
Path("ref.fa").write_text(f">chr1\n{reference}\n")
pysam.faidx("ref.fa")
Path("loci.tsv").write_text(
    "locus_id\tchr\tstart\tend\tcoordinate_system\treference_build\tprimary_motif\n"
    "L1\tchr1\t200\t236\t0-based-half-open\ttest-build\tCAG\n"
)
header = {
    "HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(reference)}],
    "RG": [{"ID": "rg1", "SM": "S1", "PL": "ONT",
            "DS": "basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.2.0"}],
}
for ext, mode in (("bam", "wb"), ("cram", "wc")):
    with pysam.AlignmentFile(f"reads.{ext}", mode, header=header, reference_filename="ref.fa") as out:
        for i in range(20):
            hp = i % 2 + 1
            read = pysam.AlignedSegment(out.header)
            read.query_name = f"read{i}"
            read.query_sequence = left + "TAC" + "CAG" * (10 if hp == 1 else 15) + "CCG" + right
            read.query_qualities = [40] * len(read.query_sequence)
            read.reference_id, read.reference_start, read.mapping_quality = 0, 0, 60
            read.cigarstring = "436M" if hp == 1 else "218M15I218M"
            read.set_tags([("RG", "rg1"), ("HP", hp), ("PS", 42)])
            out.write(read)
    pysam.index(f"reads.{ext}")
PY
python3 /pipeline/bin/catalog_adapter.py --callers medaka --locus-manifest loci.tsv \
  --reference-fai ref.fa.fai --output-dir catalogs
for ext in bam cram; do
  index=bai
  if [[ "$ext" == cram ]]; then index=crai; fi
  python3 /pipeline/bin/run_caller.py --caller medaka --sample-id S1 --platform ont --sex XY \
    --alignment "reads.$ext" --alignment-index "reads.$ext.$index" --reference ref.fa \
    --catalog catalogs/medaka.bed --config /pipeline/configs/medaka.yml --output-dir "$ext" \
    || { cat "$ext/caller.stderr.log" >&2; exit 1; }
  python3 /pipeline/bin/normalize_calls.py --caller medaka --sample-id S1 --platform ont \
    --vcf "$ext/calls.vcf.gz" --source-vcf "$ext/calls.vcf.gz" --locus-manifest loci.tsv \
    --output "$ext.tsv.gz" --versions "$ext.versions.yml"
done
python3 - <<'PY'
import csv
import gzip
for ext in ("bam", "cram"):
    with gzip.open(f"{ext}.tsv.gz", "rt") as inp:
        row, = csv.DictReader(inp, delimiter="\t")
    assert (row["canonical_locus_id"], row["call_status"]) == ("L1", "valid_call"), row
    assert sorted(int(row[f"allele{i}_length_bp"]) for i in (1, 2)) == [36, 51], row
PY
BASH
