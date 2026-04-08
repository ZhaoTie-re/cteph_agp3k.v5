#!/bin/zsh
# =============================================================================
# Prepare genotype for fixed-model fitting
# =============================================================================
# Steps:
# 1) Remove samples marked SELECTED_FOR_REMOVAL=true in PI_HAT vertex table
#    (IID-based; write FID IID without header, FID defaults to IID).
# 2) Remove monomorphic variants after sample removal (plink2 --mac 1).
# 3) Use case/ctrl/all subset to estimate MAF, then split variants by threshold.
#    - Reference subset for MAF estimation:
#      case -> FAM phenotype column == 2
#      ctrl -> FAM phenotype column == 1
#      all  -> all samples
# =============================================================================
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 <popgmm_bfile_prefix> <pihat_vertex_tsv> <out_prefix> <maf_group:ctrl|case|all> <maf_threshold> <threads>" >&2
  exit 1
fi

popgmm_bfile_prefix="$1"
pihat_vertex_tsv="$2"
out_prefix="$3"
maf_group_raw="$4"
maf_threshold="$5"
threads="$6"

maf_group="$(echo "$maf_group_raw" | tr '[:upper:]' '[:lower:]')"
if [[ "$maf_group" != "ctrl" && "$maf_group" != "case" && "$maf_group" != "all" ]]; then
  echo "[ERROR] Invalid maf_group: $maf_group_raw (allowed: ctrl, case, all)" >&2
  exit 1
fi

if ! [[ "$maf_threshold" =~ '^[0-9]*\.?[0-9]+$' ]]; then
  echo "[ERROR] Invalid maf_threshold: $maf_threshold" >&2
  exit 1
fi

log_file="${out_prefix}.fixed_model_prep.log.txt"
remove_fid_iid="${out_prefix}.pihat_selected.exclude.fid_iid"
post_remove_prefix="${out_prefix}.rm_pihat"
fixed_ready_prefix="${out_prefix}.fixed_ready"
maf_group_keep="${out_prefix}.maf_group.keep.fid_iid"
maf_ref_prefix="${out_prefix}.maf_ref"
maf_ge_list="${out_prefix}.maf_ge_threshold.variants.txt"
maf_lt_list="${out_prefix}.maf_lt_threshold.variants.txt"
maf_ge_prefix="${out_prefix}.maf_ge_threshold"
maf_lt_prefix="${out_prefix}.maf_lt_threshold"

tmp_dir="tmp"
mkdir -p "$tmp_dir"
work_iid="${tmp_dir}/${out_prefix}.work.iid.txt"
selected_iid="${tmp_dir}/${out_prefix}.pihat.selected.iid.txt"

# Build IID universe from current PopGMM genotype
awk '{ print $2 }' "${popgmm_bfile_prefix}.fam" | sed '/^$/d' | sort -u > "$work_iid"

# Parse PI_HAT vertex table: keep rows with SELECTED_FOR_REMOVAL=true-like values
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
' "$pihat_vertex_tsv" | sed '/^$/d' | sort -u > "$selected_iid"

# Restrict removal IDs to samples present in this genotype
comm -12 "$work_iid" "$selected_iid" | awk '{ print $1"\t"$1 }' > "$remove_fid_iid"

# Step 1: remove selected PI_HAT samples (skip if no overlap)
if [[ -s "$remove_fid_iid" ]]; then
  plink2 \
    --bfile "$popgmm_bfile_prefix" \
    --remove "$remove_fid_iid" \
    --make-bed \
    --out "$post_remove_prefix" \
    --threads "$threads"
else
  cp "${popgmm_bfile_prefix}.bed" "${post_remove_prefix}.bed"
  cp "${popgmm_bfile_prefix}.bim" "${post_remove_prefix}.bim"
  cp "${popgmm_bfile_prefix}.fam" "${post_remove_prefix}.fam"
fi

# Step 2: remove monomorphic variants
plink2 \
  --bfile "$post_remove_prefix" \
  --mac 1 \
  --make-bed \
  --out "$fixed_ready_prefix" \
  --threads "$threads"

# Build MAF reference subset by FAM phenotype
if [[ "$maf_group" == "ctrl" ]]; then
  awk '$6 == 1 { print $1"\t"$2 }' "${fixed_ready_prefix}.fam" > "$maf_group_keep"
elif [[ "$maf_group" == "case" ]]; then
  awk '$6 == 2 { print $1"\t"$2 }' "${fixed_ready_prefix}.fam" > "$maf_group_keep"
else
  awk '{ print $1"\t"$2 }' "${fixed_ready_prefix}.fam" > "$maf_group_keep"
fi

if [[ ! -s "$maf_group_keep" ]]; then
  echo "[ERROR] No samples available in maf_group=$maf_group under ${fixed_ready_prefix}.fam" >&2
  exit 1
fi

# Step 3: estimate AF/MAF in selected group, then split variants by threshold
plink2 \
  --bfile "$fixed_ready_prefix" \
  --keep "$maf_group_keep" \
  --freq \
  --out "$maf_ref_prefix" \
  --threads "$threads"

: > "$maf_ge_list"
: > "$maf_lt_list"

awk -v thr="$maf_threshold" '
NR == 1 {
  for (i = 1; i <= NF; i++) {
    if ($i == "ID") id_col = i
    if ($i == "ALT_FREQS") af_col = i
  }
  next
}
{
  if (id_col == 0 || af_col == 0) next
  f = $af_col
  if (f == "." || f == "NA" || f == "nan") next
  split(f, arr, ",")
  af = arr[1] + 0.0
  maf = af
  if (maf > 0.5) maf = 1.0 - maf
  if (maf >= thr) {
    print $id_col >> ge_file
  } else {
    print $id_col >> lt_file
  }
}
' ge_file="$maf_ge_list" lt_file="$maf_lt_list" "${maf_ref_prefix}.afreq"

if [[ ! -s "$maf_ge_list" ]]; then
  echo "[ERROR] No variants with MAF >= $maf_threshold in maf_group=$maf_group" >&2
  exit 1
fi

plink2 \
  --bfile "$fixed_ready_prefix" \
  --extract "$maf_ge_list" \
  --make-bed \
  --out "$maf_ge_prefix" \
  --threads "$threads"

if [[ -s "$maf_lt_list" ]]; then
  plink2 \
    --bfile "$fixed_ready_prefix" \
    --extract "$maf_lt_list" \
    --make-bed \
    --out "$maf_lt_prefix" \
    --threads "$threads"
fi

before_n=$(wc -l < "${popgmm_bfile_prefix}.fam")
before_v=$(wc -l < "${popgmm_bfile_prefix}.bim")
after_rm_n=$(wc -l < "${post_remove_prefix}.fam")
after_rm_v=$(wc -l < "${post_remove_prefix}.bim")
after_fix_n=$(wc -l < "${fixed_ready_prefix}.fam")
after_fix_v=$(wc -l < "${fixed_ready_prefix}.bim")
ge_n=$(wc -l < "$maf_ge_list")
lt_n=$(wc -l < "$maf_lt_list")

{
  echo "[$(date)] Fixed-model genotype preparation summary"
  echo "INPUT_BFILE_PREFIX: $popgmm_bfile_prefix"
  echo "PIHAT_VERTEX_TSV: $pihat_vertex_tsv"
  echo "MAF_GROUP: $maf_group"
  echo "MAF_THRESHOLD: $maf_threshold"
  echo ""
  echo "Generated files:"
  echo "  remove list           : $remove_fid_iid"
  echo "  fixed-ready bfile     : ${fixed_ready_prefix}.{bed,bim,fam}"
  echo "  maf group keep list   : $maf_group_keep"
  echo "  af/maf reference      : ${maf_ref_prefix}.afreq"
  echo "  MAF >= threshold list : $maf_ge_list"
  echo "  MAF <  threshold list : $maf_lt_list"
  echo "  MAF >= threshold bfile: ${maf_ge_prefix}.{bed,bim,fam}"
  if [[ -s "$maf_lt_list" ]]; then
    echo "  MAF <  threshold bfile: ${maf_lt_prefix}.{bed,bim,fam}"
  else
    echo "  MAF <  threshold bfile: not generated (0 variants)"
  fi
  echo ""
  echo "Counts (samples / variants):"
  echo "  Input genotype            : ${before_n} / ${before_v}"
  echo "  After PIHAT sample rm     : ${after_rm_n} / ${after_rm_v}"
  echo "  After monomorphic rm      : ${after_fix_n} / ${after_fix_v}"
  echo ""
  echo "MAF split variant counts:"
  echo "  MAF >= threshold          : ${ge_n}"
  echo "  MAF <  threshold          : ${lt_n}"
} > "$log_file"

echo "[OK] Fixed-model genotype prepared: $out_prefix"
