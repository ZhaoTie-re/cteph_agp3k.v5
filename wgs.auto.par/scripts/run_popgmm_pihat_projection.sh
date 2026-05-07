#!/bin/zsh
# =============================================================================
# PopGMM-PIHAT PCA Base and Projection Pipeline
# =============================================================================
# Purpose:
#   Build a PCA base after excluding samples in the overlap between
#   PopGMM samples and PI_HAT SELECTED_FOR_REMOVAL samples, then project
#   PopGMM samples onto that PCA base.
#
# Inputs:
#   1) popgmm_fid_iid         : PopGMM sample list (FID IID, no header)
#   2) pihat_vertex_tsv       : PI_HAT vertex-cover TSV with IID and
#                               SELECTED_FOR_REMOVAL columns
#   3) popgmm_bfile_prefix    : PLINK prefix (bed/bim/fam)
#   4) high_ld_regions        : BED/range file for high-LD exclusion
#   5) threads                : Number of computation threads
#   6) out_prefix             : Output prefix
#   7) cleanup_tmp (optional) : true/false, default false
#
# Outputs:
#   - overlap exclude list (FID IID)
#   - base genotype / prune list / PCA artifacts (eigenval/eigenvec/acount)
#   - projected scores for PopGMM samples (.sscore/.sscore.vars)
#   - done marker file
#
# Workflow:
#   Step 1) Build overlap exclude list from PopGMM and PI_HAT removals
#   Step 2) Build PCA base on non-overlap samples
#   Step 3) Project PopGMM samples using base allele weights/frequencies
# =============================================================================
set -euo pipefail

if [[ $# -lt 6 || $# -gt 7 ]]; then
  echo "Usage: $0 <popgmm_fid_iid> <pihat_vertex_tsv> <popgmm_bfile_prefix> <high_ld_regions> <threads> <out_prefix> [cleanup_tmp=true|false]" >&2
  exit 1
fi

popgmm_list_fid_iid="$1"
pihat_vertex_cover_tsv="$2"
popgmm_genotype_prefix="$3"
high_ld_regions_bed="$4"
num_threads="$5"
output_prefix="$6"
cleanup_tmp="${7:-false}"

tmp_dir="tmp"
mkdir -p "${tmp_dir}"

popgmm_iid_list="${tmp_dir}/${output_prefix}.popgmm.iids.txt"
pihat_selected_removal_iids="${tmp_dir}/${output_prefix}.pihat_selected_for_removal.iids.txt"
intersection_iid_list="${output_prefix}.pihat_popgmm_intersection.iid"
overlap_exclude_fid_iid="${output_prefix}.pihat_popgmm_overlap.exclude.fid_iid"

# -----------------------------------------------------------------------------
# Step 1. Build intersection sample list
# - Input A: PopGMM list (FID IID without header; IID-only also tolerated)
# - Input B: PI_HAT vertex-cover table rows with SELECTED_FOR_REMOVAL=True
# - Output : FID/IID exclude file where FID == IID
# -----------------------------------------------------------------------------
awk '{ if (NF >= 2) print $2; else if (NF >= 1) print $1; }' "${popgmm_list_fid_iid}" | sed '/^$/d' | sort -u > "${popgmm_iid_list}"

awk -F'\t' '
BEGIN { iid_col = 0; sel_col = 0 }
NR == 1 {
  for (i = 1; i <= NF; i++) {
    if ($i == "IID") iid_col = i
    if ($i == "SELECTED_FOR_REMOVAL") sel_col = i
  }
  next
}
{
  if (iid_col == 0 || sel_col == 0) next
  v = tolower($sel_col)
  if (v == "true" || v == "1" || v == "yes") print $iid_col
}
' "${pihat_vertex_cover_tsv}" | sed '/^$/d' | sort -u > "${pihat_selected_removal_iids}"

comm -12 "${popgmm_iid_list}" "${pihat_selected_removal_iids}" > "${intersection_iid_list}"
awk '{ print $1"\t"$1 }' "${intersection_iid_list}" > "${overlap_exclude_fid_iid}"

if [[ ! -s "${overlap_exclude_fid_iid}" ]]; then
  echo "[ERROR] Intersection between PopGMM and PI_HAT selected-for-removal IIDs is empty: ${overlap_exclude_fid_iid}" >&2
  exit 1
fi

# -----------------------------------------------------------------------------
# Step 2. Build PCA base from non-intersection samples
# - Remove intersection samples from PopGMM genotype
# - LD pruning on autosomal common SNPs outside high-LD regions
# - Compute PCA with allele weights and allele counts
# -----------------------------------------------------------------------------
base_samples_prefix="${output_prefix}.base_no_intersection"
base_prune_prefix="${output_prefix}.base_no_intersection.prune"
base_pca_prefix="${output_prefix}.base_no_intersection.pca"

plink2 \
  --bfile "${popgmm_genotype_prefix}" \
  --remove "${overlap_exclude_fid_iid}" \
  --make-bed \
  --out "${base_samples_prefix}" \
  --threads "${num_threads}"

plink2 \
  --bfile "${base_samples_prefix}" \
  --autosome \
  --snps-only just-acgt \
  --maf 0.01 \
  --exclude range "${high_ld_regions_bed}" \
  --indep-pairwise 50 5 0.2 \
  --out "${base_prune_prefix}" \
  --threads "${num_threads}"

plink2 \
  --bfile "${base_samples_prefix}" \
  --extract "${base_prune_prefix}.prune.in" \
  --freq counts \
  --pca 20 allele-wts approx \
  --out "${base_pca_prefix}" \
  --threads "${num_threads}"

# -----------------------------------------------------------------------------
# Step 3. Project all PopGMM samples onto the PCA base
# - Restrict all samples to the same pruned SNP set from Step 2
# - Project using base allele frequencies and eigenvector allele weights
# -----------------------------------------------------------------------------
projection_input_prefix="${output_prefix}.all_samples_for_projection"
projection_score_prefix="${output_prefix}.all_samples_projection"

# Keep this genotype as a temporary intermediate by default under tmp/.
projection_input_prefix="${tmp_dir}/${projection_input_prefix}"

plink2 \
  --bfile "${popgmm_genotype_prefix}" \
  --extract "${base_prune_prefix}.prune.in" \
  --make-bed \
  --out "${projection_input_prefix}" \
  --threads "${num_threads}"

plink2 \
  --bfile "${projection_input_prefix}" \
  --read-freq "${base_pca_prefix}.acount" \
  --score "${base_pca_prefix}.eigenvec.allele" 2 6 header-read no-mean-imputation variance-standardize list-variants \
  --score-col-nums 7-26 \
  --out "${projection_score_prefix}" \
  --threads "${num_threads}"

if [[ "${cleanup_tmp}" == "true" ]]; then
  rm -rf "${tmp_dir}"
fi

echo "[OK] PopGMM-PIHAT projection completed: ${output_prefix}" > "${output_prefix}.done.txt"
