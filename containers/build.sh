#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root_dir"

docker build --pull --file containers/Dockerfile.normalizer --tag local/strpadre-normalizer:1.0.0 .
docker build --file containers/Dockerfile.longtr --tag local/strpadre-longtr:1.2 .
docker build --file containers/Dockerfile.atarva --tag local/strpadre-atarva:0.7.1 .
docker build --file containers/Dockerfile.strdust --tag local/strpadre-strdust:0.20.0 .

if [[ -f containers/vendor/trgt-5.1.0 ]]; then
  docker build --file containers/Dockerfile.trgt --tag local/strpadre-trgt:5.1.0 .
else
  printf '%s\n' 'TRGT image not built: place an authorised 5.1.0 binary at containers/vendor/trgt-5.1.0 first.' >&2
fi
