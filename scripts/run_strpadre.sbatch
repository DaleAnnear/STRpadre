#!/usr/bin/env bash
#SBATCH --job-name=strpadre
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=strpadre-%j.out
#SBATCH --error=strpadre-%j.err

# Submit from the repository root:
#   PARAMS_FILE=/path/to/params.yml sbatch scripts/run_strpadre.sbatch
#
# For a HiFi run with a locally built, authorised TRGT image:
#   CALLERS=trgt,longtr,atarva,strdust \
#     PARAMS_FILE=/path/to/params.yml sbatch scripts/run_strpadre.sbatch
#
# Optional site-specific scheduler settings:
##SBATCH --partition=CHANGE_ME
##SBATCH --account=CHANGE_ME

set -euo pipefail

# Uncomment and adapt these lines if your site provides Nextflow and Docker as modules.
# module load nextflow
# module load docker

PROJECT_DIR="${PROJECT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}"
PARAMS_FILE="${PARAMS_FILE:-${PROJECT_DIR}/params.yml}"
PROFILE="${PROFILE:-docker}"
CALLERS="${CALLERS:-longtr,atarva,strdust}"
MAX_FORKS="${MAX_FORKS:-1}"

# Keep Nextflow metadata and work files on shared storage by default so -resume works
# across SLURM jobs. Override NXF_WORK with a shared high-performance filesystem if
# your site provides one.
export NXF_HOME="${NXF_HOME:-${PROJECT_DIR}/.nextflow}"
export NXF_WORK="${NXF_WORK:-${PROJECT_DIR}/work}"

if [[ ! "${MAX_FORKS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "MAX_FORKS must be a positive integer; received: ${MAX_FORKS}" >&2
    exit 2
fi

for executable in nextflow docker; do
    if ! command -v "${executable}" >/dev/null 2>&1; then
        echo "Required executable not found: ${executable}" >&2
        exit 127
    fi
done

if [[ ! -f "${PARAMS_FILE}" ]]; then
    echo "Parameter file not found: ${PARAMS_FILE}" >&2
    exit 2
fi

if ! docker info >/dev/null 2>&1; then
    echo "Docker is unavailable to this SLURM job. Load/configure Docker or use a site-specific container profile." >&2
    exit 1
fi

mkdir -p "${NXF_HOME}" "${NXF_WORK}"
cd "${PROJECT_DIR}"

# The checked-in Docker profile uses UID/GID 1000. Override it with the submitting
# user's identity so mounted results and work files remain writable. maxForks=1 keeps
# the local Nextflow executor within this single SLURM allocation by default.
RUNTIME_CONFIG="$(mktemp "${NXF_HOME}/strpadre-${SLURM_JOB_ID:-local}.XXXXXX.config")"
trap 'rm -f "${RUNTIME_CONFIG}"' EXIT
cat > "${RUNTIME_CONFIG}" <<EOF
process {
    maxForks = ${MAX_FORKS}
}

docker.runOptions = '--user $(id -u):$(id -g)'
EOF

echo "SLURM job: ${SLURM_JOB_ID:-interactive}"
echo "Project: ${PROJECT_DIR}"
echo "Parameters: ${PARAMS_FILE}"
echo "Profile: ${PROFILE}"
echo "Callers: ${CALLERS}"
echo "Nextflow work directory: ${NXF_WORK}"

nextflow run "${PROJECT_DIR}/main.nf" \
    -profile "${PROFILE}" \
    -params-file "${PARAMS_FILE}" \
    -c "${RUNTIME_CONFIG}" \
    --callers "${CALLERS}" \
    -resume \
    "$@"
