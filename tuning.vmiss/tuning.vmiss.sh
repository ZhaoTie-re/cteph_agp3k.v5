#!/bin/zsh
#SBATCH --job-name=tuning_vmiss_pipeline
#SBATCH --output=tuning_vmiss_%j.out
#SBATCH --error=tuning_vmiss_%j.err
#SBATCH -p gr10478b
#SBATCH -t 168:0:0
#SBATCH --rsc p=1:t=6:c=6

# ============================================================================
# VMISS Threshold Tuning Pipeline - SLURM Batch Script
# ============================================================================
# Purpose: Identify optimal variant missingness (VMISS) thresholds by
#          sequencing depth and/or allele frequency stratification.
#
# Requirements:
#   - Conda environment: cteph_geno_pro
#   - Tools: python3, plink2
#   - Input: PLINK bed/bim/fam files + metadata xlsx
#
# Output:
#   - variant_metrics_compiled.txt        : VMISS/MAF metrics
#   - vmiss_threshold_analysis.json       : Analysis results
#   - vmiss_threshold_recommendations.tsv : Threshold recommendations
#   - vmiss_*_analysis_modeX.png          : Distribution plots
#   - run_manifest.json                   : Execution metadata
#   - run_vmiss_tuning.log                : Pipeline log
# ============================================================================
set -euo pipefail

log_info()  { echo "[INFO] $*"; }
log_warn()  { echo "[WARN] $*"; }
log_error() { echo "[ERROR] $*" >&2; }
die()       { log_error "$*"; exit 1; }

require_cmd() {
    local cmd="$1"
    command -v "$cmd" >/dev/null 2>&1 || die "Required command not found: ${cmd}"
}

is_true() {
    case "${1:l}" in
        true|1|yes|y) return 0 ;;
        *) return 1 ;;
    esac
}

# ============================================================================
# RUNTIME INITIALIZATION
# ============================================================================
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    # Running under SLURM
    WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
    log_info "Detected SLURM job. Using submission directory: ${WORKDIR}"
else
    # Running interactively
    if [[ -n "${ZSH_VERSION:-}" ]]; then
        SCRIPT_SOURCE="${(%):-%x}"
    else
        SCRIPT_SOURCE="${BASH_SOURCE[0]}"
    fi
    SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_SOURCE}")" && pwd)"
    WORKDIR="${SCRIPT_DIR}"
fi
cd "${WORKDIR}" || die "Could not change to ${WORKDIR}"
log_info "Working directory: ${WORKDIR}"

PYTHON_SCRIPT="${WORKDIR}/scripts/run_vmiss_tuning.py"
[[ -f "${PYTHON_SCRIPT}" ]] || die "Pipeline entry script not found: ${PYTHON_SCRIPT}"

# Activate Conda environment
log_info "Loading conda environment..."
require_cmd conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cteph_geno_pro || die "Failed to activate conda env: cteph_geno_pro"

# Validate dependencies
export PATH="/home/b/b37974/plink2_alpha6:${PATH}"
require_cmd plink2
require_cmd python3
log_info "Dependencies verified"

# ============================================================================
# PARAMETERS
# ============================================================================
# User-configurable parameters (environment-variable override supported)
# Usage example:
#   ID_COL="ID JHRPv5" THREADS=8 SKIP_CALC=false sbatch tuning.vmiss.sh

# [Execution control]
# SKIP_CALC=true  => skip metric calculation and analyze existing metrics file
# SKIP_CALC=false => run full pipeline (calculation + analysis)
SKIP_CALC="${SKIP_CALC:-true}"

# [I/O Paths]
# BED_PREFIX: PLINK prefix (without .bed/.bim/.fam)
# INFO_PATH : Metadata file path (.xlsx/.csv/.tsv/.txt)
# OUT_DIR   : Output directory
BED_PREFIX="${BED_PREFIX:-/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/wgs.auto.par/results/07_sample_qc/run_qc/cteph_agp3k_v5_wgs_merged.sample_qc}"
INFO_PATH="${INFO_PATH:-/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/info/ph_agp3k_combined_final_20251112.xlsx}"
OUT_DIR="${OUT_DIR:-/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/tuning.vmiss/results}"

# [Metadata schema]
# ID_COL    : Sample ID column in metadata; matched against FAM IID
# DP_COL    : Sequencing depth column (used for 15X/30X grouping)
# PHENO_COL : Phenotype/group column
# CASE_VALUE/CTRL_VALUE: Values in PHENO_COL defining case/control groups
ID_COL="${ID_COL:-ID JHRPv5}"
DP_COL="${DP_COL:-Target DP}"
PHENO_COL="${PHENO_COL:-Outcome}"
CASE_VALUE="${CASE_VALUE:-PH}"
CTRL_VALUE="${CTRL_VALUE:-AGP3K}"

# [Compute]
# THREADS: Threads passed to Python pipeline and downstream plink2 calls.
#          Defaults to SLURM_CPUS_PER_TASK when available, otherwise 16.
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"

# [MAF calculation strategy]
# MAF_GROUP: ctrl (control-only), case (case-only), all (pooled)
MAF_GROUP="${MAF_GROUP:-ctrl}"

# [Analysis mode]
# 1: Unified
# 2: Depth-Stratified (recommended)
# 3: MAF-Stratified
# 4: Combined (Depth × MAF)
ANALYSIS_MODE="${ANALYSIS_MODE:-2}"

# [CDF resolution]
# Smaller step => finer resolution but potentially slower computation
CDF_STEP="${CDF_STEP:-0.001}"
CDF_STEP_15X="${CDF_STEP_15X:-0.01}"
CDF_STEP_30X="${CDF_STEP_30X:-0.01}"

# [Kneedle algorithm weights]
# KNEE_WEIGHT_X/KNEE_WEIGHT_Y tune sensitivity along x/y axes
# v4-equivalent setting: KNEE_WEIGHT_X=4 and KNEE_WEIGHT_Y=1
# Current run: use slightly looser threshold to retain more variants (especially around RNF213)
# Practical knob: decrease KNEE_WEIGHT_X (e.g., 2) to make inclusion less strict
KNEE_WEIGHT_X="${KNEE_WEIGHT_X:-2}"
KNEE_WEIGHT_Y="${KNEE_WEIGHT_Y:-1}"

# [Output threshold rounding]
THRESHOLD_DECIMALS="${THRESHOLD_DECIMALS:-2}"

SKIP_CALC_FLAG=""
if is_true "${SKIP_CALC}"; then
    SKIP_CALC_FLAG="--skip-calc"
fi

mkdir -p "${OUT_DIR}"

echo ""
echo "========== PREFLIGHT VALIDATION =========="

# Check input files
for ext in bed bim fam; do
    [[ -f "${BED_PREFIX}.${ext}" ]] || die "Missing input file: ${BED_PREFIX}.${ext}"
done
log_info "PLINK files verified"

[[ -f "${INFO_PATH}" ]] || die "Metadata file not found: ${INFO_PATH}"
log_info "Metadata file verified"

[[ "${THREADS}" =~ ^[0-9]+$ && "${THREADS}" -gt 0 ]] || die "Invalid THREADS value: ${THREADS}"
[[ "${ANALYSIS_MODE}" =~ ^[1-4]$ ]] || die "Invalid ANALYSIS_MODE (1-4): ${ANALYSIS_MODE}"
log_info "Parameters validated"

if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]] && [[ "${THREADS}" -gt "${SLURM_CPUS_PER_TASK}" ]]; then
    log_warn "THREADS=${THREADS} exceeds SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}; Python will clamp to available CPUs"
fi

echo ""
echo "========== CONFIGURATION =========="
echo "Dataset:   ${BED_PREFIX}"
echo "Metadata:  ${INFO_PATH}"
echo "Output:    ${OUT_DIR}"
echo "Threads:   ${THREADS}"
echo "ID Col:    ${ID_COL}"
echo "Mode:      ${ANALYSIS_MODE} | MAF: ${MAF_GROUP}"
echo ""

START_TS=$(date +%s)
echo "========== EXECUTION =========="
echo ""

cmd=(
    python3 "${PYTHON_SCRIPT}"
    --bfile "${BED_PREFIX}"
    --info "${INFO_PATH}"
    --out-dir "${OUT_DIR}"
    --id-col "${ID_COL}"
    --dp-col "${DP_COL}"
    --pheno-col "${PHENO_COL}"
    --case-value "${CASE_VALUE}"
    --ctrl-value "${CTRL_VALUE}"
    --maf-group "${MAF_GROUP}"
    --analysis-mode "${ANALYSIS_MODE}"
    --cdf-step "${CDF_STEP}"
    --cdf-step-15x "${CDF_STEP_15X}"
    --cdf-step-30x "${CDF_STEP_30X}"
    --knee-weight-x "${KNEE_WEIGHT_X}"
    --knee-weight-y "${KNEE_WEIGHT_Y}"
    --threshold-decimals "${THRESHOLD_DECIMALS}"
    --threads "${THREADS}"
    --log-file "${OUT_DIR}/run_vmiss_tuning.log"
)

if [[ -n "${SKIP_CALC_FLAG}" ]]; then
    cmd+=("${SKIP_CALC_FLAG}")
fi

set +e
"${cmd[@]}"
PIPELINE_EXIT_CODE=$?
set -e

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))

echo ""
echo "========== PIPELINE SUMMARY =========="
if [[ ${PIPELINE_EXIT_CODE} -eq 0 ]]; then
    echo "Status:         SUCCESS"
    echo "Execution Time: ${ELAPSED}s"
    echo "Output Dir:     ${OUT_DIR}"
    echo ""
    echo "Generated Files:"
    echo "  - variant_metrics_compiled.txt"
    echo "  - vmiss_threshold_recommendations.tsv"
    echo "  - vmiss_threshold_analysis.json"
    echo "  - vmiss_*_analysis_mode${ANALYSIS_MODE}.png"
    echo "  - run_manifest.json"
    echo "  - run_vmiss_tuning.log"
else
    echo "Status:         FAILED"
    echo "Exit Code:      ${PIPELINE_EXIT_CODE}"
    echo "Execution Time: ${ELAPSED}s"
    echo ""
    echo "Troubleshooting:"
    echo "  1. Check ${OUT_DIR}/run_vmiss_tuning.log for error messages"
    echo "  2. Verify all input files exist and are readable"
    echo "  3. Review run_manifest.json for parameter configuration"
    echo "  4. Ensure conda environment 'cteph_geno_pro' is properly installed"
    exit ${PIPELINE_EXIT_CODE}
fi
echo ""
