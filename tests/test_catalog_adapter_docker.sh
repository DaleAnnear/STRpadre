#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root_dir"
runner_uid="$(id -u):$(id -g)"
docker build --file containers/Dockerfile.normalizer --tag local/strpadre-normalizer:1.0.0 .

tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT
cp examples/loci.tsv "$tmp_dir/loci.tsv"
cp examples/reference/tiny.fa.fai "$tmp_dir/reference.fa.fai"
docker run --rm --user "$runner_uid" -v "$root_dir:/workspace:ro" -v "$tmp_dir:/run" -w /run local/strpadre-normalizer:1.0.0 \
  python /workspace/bin/catalog_adapter.py --locus-manifest loci.tsv --reference-fai reference.fa.fai \
    --callers trgt,longtr,atarva,strdust --output-dir catalogs
for output in trgt.bed longtr.bed atarva.bed.gz atarva.bed.gz.tbi strdust.bed catalog_locus_mapping.tsv catalog_adaptation_report.json; do
  test -s "$tmp_dir/catalogs/$output"
done
