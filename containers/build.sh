#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root_dir"

usage() {
  cat <<'EOF'
Usage: bash containers/build.sh [--local] [IMAGE ...]

Build or acquire STRpadre container images.

With no IMAGE arguments, all images are selected. By default, published images
are pulled from Docker Hub and built locally only when the pull fails. --local
forces local builds for the selected images.

Images: normalizer, longtr, atarva, strdust, trgt, all
EOF
}

force_local=false
selected=()

add_image() {
  local image=$1
  local existing
  for existing in "${selected[@]}"; do
    [[ "$existing" == "$image" ]] && return
  done
  selected+=("$image")
}

while (( $# )); do
  case "$1" in
    --local)
      force_local=true
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    all)
      add_image normalizer
      add_image longtr
      add_image atarva
      add_image strdust
      add_image trgt
      ;;
    normalizer|longtr|atarva|strdust|trgt)
      add_image "$1"
      ;;
    *)
      printf 'Unknown option or image: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
  shift
done

if (( ${#selected[@]} == 0 )); then
  selected=(normalizer longtr atarva strdust trgt)
fi

normalizer_local=local/strpadre-normalizer:1.0.0
normalizer_hub=daleannear/strpadre:normalizer-1.0.0

pull_or_build() {
  local local_image=$1
  local hub_image=$2
  local dockerfile=$3
  local requires_normalizer=${4:-false}

  if [[ "$force_local" == false ]] && docker pull "$hub_image"; then
    docker tag "$hub_image" "$local_image"
    return
  fi

  if [[ "$force_local" == true ]]; then
    printf 'Building %s locally.\n' "$local_image" >&2
  else
    printf 'Docker Hub image unavailable; building %s locally.\n' "$local_image" >&2
  fi
  if [[ "$requires_normalizer" == true ]]; then
    ensure_local_normalizer
  fi
  docker build --pull --file "$dockerfile" --tag "$local_image" .
}

ensure_local_normalizer() {
  if ! docker image inspect "$normalizer_local" >/dev/null 2>&1; then
    pull_or_build "$normalizer_local" "$normalizer_hub" containers/Dockerfile.normalizer
  fi
}

for image in "${selected[@]}"; do
  case "$image" in
    normalizer)
      pull_or_build "$normalizer_local" "$normalizer_hub" containers/Dockerfile.normalizer
      ;;
    longtr)
      pull_or_build local/strpadre-longtr:1.2 daleannear/strpadre:longtr-1.2 containers/Dockerfile.longtr true
      ;;
    atarva)
      pull_or_build local/strpadre-atarva:0.7.1 daleannear/strpadre:atarva-0.7.1 containers/Dockerfile.atarva true
      ;;
    strdust)
      pull_or_build local/strpadre-strdust:0.20.0 daleannear/strpadre:strdust-0.20.0 containers/Dockerfile.strdust true
      ;;
    trgt)
      if [[ -f containers/vendor/trgt-5.1.0 ]]; then
        ensure_local_normalizer
        docker build --file containers/Dockerfile.trgt --tag local/strpadre-trgt:5.1.0 .
      else
        printf '%s\n' 'TRGT image not built: place an authorised 5.1.0 binary at containers/vendor/trgt-5.1.0 first.' >&2
      fi
      ;;
  esac
done
