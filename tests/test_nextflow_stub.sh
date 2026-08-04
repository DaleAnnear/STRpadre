#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root_dir"
runner_uid="$(id -u):$(id -g)"

nextflow_cmd=(docker run --rm -e HOME=/tmp -e NXF_HOME=/tmp/nxf -e NXF_TEMP=/tmp --user "$runner_uid" -v "$root_dir:/workspace" -w /workspace nextflow/nextflow:24.10.4 nextflow)
for callers in trgt longtr atarva strdust trgt,longtr,atarva,strdust; do
  "${nextflow_cmd[@]}" run main.nf -c conf/test.config -profile test -stub-run -params-file tests/data/stub.params.yml --callers "$callers" --outdir /tmp/strpadre-test-output -work-dir "/tmp/strpadre-work-${callers//,/_}"
done

if "${nextflow_cmd[@]}" run main.nf -c conf/test.config -profile test -stub-run -params-file tests/data/stub.params.yml --callers badcaller --outdir /tmp/strpadre-test-output -work-dir /tmp/strpadre-work-invalid; then
  printf '%s\n' 'invalid caller unexpectedly launched' >&2
  exit 1
fi
