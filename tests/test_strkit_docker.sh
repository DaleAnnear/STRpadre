#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
docker run --rm -i -v "$root_dir:/pipeline:ro" -w /tmp local/strpadre-strkit:0.25.0 bash -eu <<'BASH'
python - <<'PY'
import random
from pathlib import Path
import pysam

rng = random.Random(12)
flank = lambda: "".join(rng.choices("ACGT", k=200))
seq = flank() + "CAG" * 10 + flank()
Path("ref.fa").write_text(f">chr1\n{seq}\n")
pysam.faidx("ref.fa")
Path("loci.tsv").write_text(
    "locus_id\tchr\tstart\tend\tcoordinate_system\treference_build\tprimary_motif\n"
    "L1\tchr1\t200\t230\t0-based-half-open\ttest-build\tCAG\n"
)
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(seq)}]}
for ext, mode in [("bam", "wb"), ("cram", "wc")]:
    name = f"reads.{ext}"
    with pysam.AlignmentFile(name, mode, header=header, reference_filename="ref.fa") as output:
        for i in range(16):
            read = pysam.AlignedSegment(output.header)
            read.query_name = f"read{i}"
            read.query_sequence = seq
            read.query_qualities = [40] * len(seq)
            read.reference_id, read.reference_start, read.mapping_quality = 0, 0, 60
            read.cigarstring = f"{len(seq)}M"
            output.write(read)
    pysam.index(name)
PY
python /pipeline/bin/catalog_adapter.py --callers strkit --locus-manifest loci.tsv \
  --reference-fai ref.fa.fai --output-dir catalogs
for platform in hifi ont; do
  alignment=reads.bam
  index=reads.bam.bai
  if [[ "$platform" == ont ]]; then
    alignment=reads.cram
    index=reads.cram.crai
  fi
  python /pipeline/bin/run_caller.py --caller strkit --sample-id S1 --platform "$platform" \
    --alignment "$alignment" --alignment-index "$index" --reference ref.fa \
    --catalog catalogs/strkit.bed --config /pipeline/configs/strkit.yml --output-dir "$platform"
  python /pipeline/bin/normalize_calls.py --caller strkit --sample-id S1 --platform "$platform" \
    --vcf "$platform/calls.vcf.gz" --source-vcf "$platform/calls.vcf.gz" --locus-manifest loci.tsv \
    --output "$platform.tsv.gz" --versions "$platform.versions.yml"
done
python - <<'PY'
import csv
import gzip
for platform in ("hifi", "ont"):
    with gzip.open(f"{platform}.tsv.gz", "rt") as f:
        row, = csv.DictReader(f, delimiter="\t")
    assert (row["canonical_locus_id"], row["call_status"]) == ("L1", "valid_call"), row
    assert (row["allele1_length_bp"], row["allele2_length_bp"]) == ("30", "30"), row
    assert (row["allele1_copy_number"], row["allele2_copy_number"]) == ("10", "10"), row
PY
BASH
