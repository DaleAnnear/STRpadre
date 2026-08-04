#!/usr/bin/env bash
set -euo pipefail

docker run --rm local/strpadre-normalizer:1.0.0 python --version
docker run --rm local/strpadre-normalizer:1.0.0 bcftools --version
docker run --rm local/strpadre-longtr:1.2 LongTR --help
docker run --rm local/strpadre-atarva:0.7.1 atarva --version
docker run --rm local/strpadre-strdust:0.20.0 STRdust --version
if docker image inspect local/strpadre-trgt:5.1.0 >/dev/null 2>&1; then
  docker run --rm local/strpadre-trgt:5.1.0 trgt --version
fi
