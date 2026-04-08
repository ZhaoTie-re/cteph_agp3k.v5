#!/bin/zsh
# =============================================================================
# Prepare genotype for random-model fitting
# =============================================================================
# Steps:
# 1) Extract variants with MAF >= threshold (from fixed-model variant list)
#    from PopGMM-subset genotype for random-effects model.
# =============================================================================
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <popgmm_bfile_prefix> <maf_ge_variants_list> <out_prefix> <threads>" >&2
  exit 1
fi

popgmm_bfile_prefix="$1"
maf_ge_variants_list="$2"
out_prefix="$3"
threads="$4"

log_file="${out_prefix}.random_model_prep.log.txt"

if [[ ! -f "${popgmm_bfile_prefix}.bed" ]]; then
  echo "[ERROR] Input bfile not found: ${popgmm_bfile_prefix}.bed" >&2
  exit 1
fi

if [[ ! -f "${maf_ge_variants_list}" ]]; then
  echo "[ERROR] MAF >= threshold variants list not found: ${maf_ge_variants_list}" >&2
  exit 1
fi

# Extract variants with MAF >= threshold for random model
plink2 \
  --bfile "${popgmm_bfile_prefix}" \
  --extract "${maf_ge_variants_list}" \
  --make-bed \
  --out "${out_prefix}" \
  --threads "${threads}"

# Generate summary log
pop_n=$(wc -l < "${popgmm_bfile_prefix}.fam")
pop_v=$(wc -l < "${popgmm_bfile_prefix}.bim")
random_n=$(wc -l < "${out_prefix}.fam")
random_v=$(wc -l < "${out_prefix}.bim")
extracted_v=$(wc -l < "${maf_ge_variants_list}")

{
  echo "[$(date)] Random model genotype preparation summary"
  echo "INPUT_BFILE_PREFIX: $popgmm_bfile_prefix"
  echo "MAF_GE_THRESHOLD_VARIANTS_LIST: $maf_ge_variants_list"
  echo ""
  echo "Counts (samples / variants):"
  echo "  PopGMM input genotype         : ${pop_n} / ${pop_v}"
  echo "  Variants to extract (MAF>=)  : ${extracted_v}"
  echo "  Random model subset           : ${random_n} / ${random_v}"
  echo ""
  echo "Summary:"
  echo "  All samples retained: ${random_n}"
  echo "  Only MAF >= threshold variants used: ${random_v}"
} > "$log_file"

echo "[OK] Random model genotype prepared: $out_prefix"
