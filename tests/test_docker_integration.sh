#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root_dir"
runner_uid="$(id -u):$(id -g)"
docker build --file containers/Dockerfile.normalizer --tag local/strpadre-normalizer:1.0.0 .

tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT
cp tests/fixtures/manifest.tsv "$tmp_dir/loci.tsv"
printf '%s\n' 'sample_id,platform,alignment,alignment_index' 'S1,hifi,a.bam,a.bai' > "$tmp_dir/samples.csv"
for caller in trgt longtr atarva strdust; do
  docker run --rm --user "$runner_uid" -v "$root_dir:/workspace:ro" -v "$tmp_dir:/run" -w /run local/strpadre-normalizer:1.0.0 \
    python /workspace/bin/normalize_calls.py --caller "$caller" --sample-id S1 --platform hifi \
      --vcf "/workspace/tests/fixtures/$caller.vcf" --locus-manifest loci.tsv --source-vcf "native/$caller.vcf" \
      --output "$caller.tsv.gz" --versions "$caller.versions.yml"
done
docker run --rm --user "$runner_uid" -v "$root_dir:/workspace:ro" -v "$tmp_dir:/run" -w /run local/strpadre-normalizer:1.0.0 \
  python /workspace/bin/build_consensus.py --normalized-files trgt.tsv.gz longtr.tsv.gz atarva.tsv.gz strdust.tsv.gz \
    --locus-manifest loci.tsv --samplesheet samples.csv --config /workspace/configs/consensus.yml \
    --callers trgt,longtr,atarva,strdust --output-dir consensus
test -s "$tmp_dir/consensus/consensus.tsv.gz"
test -s "$tmp_dir/consensus/consensus.vcf"
