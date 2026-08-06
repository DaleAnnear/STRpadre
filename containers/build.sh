#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root_dir"

pull_or_build() {
  local local_image=$1
  local hub_image=$2
  local dockerfile=$3

  if docker pull "$hub_image"; then
    docker tag "$hub_image" "$local_image"
    return
  fi

  printf 'Docker Hub image unavailable; building %s locally.\n' "$local_image" >&2
  docker build --pull --file "$dockerfile" --tag "$local_image" .
}

pull_or_build local/strpadre-normalizer:1.0.0 daleannear/strpadre:normalizer-1.0.0 containers/Dockerfile.normalizer
pull_or_build local/strpadre-longtr:1.2 daleannear/strpadre:longtr-1.2 containers/Dockerfile.longtr
pull_or_build local/strpadre-atarva:0.7.1 daleannear/strpadre:atarva-0.7.1 containers/Dockerfile.atarva
pull_or_build local/strpadre-strdust:0.20.0 daleannear/strpadre:strdust-0.20.0 containers/Dockerfile.strdust

if [[ -f containers/vendor/trgt-5.1.0 ]]; then
  docker build --file containers/Dockerfile.trgt --tag local/strpadre-trgt:5.1.0 .
else
  printf '%s\n' 'TRGT image not built: place an authorised 5.1.0 binary at containers/vendor/trgt-5.1.0 first.' >&2
fi
