#!/bin/zsh
# =============================================================================
# BBJ (Japanese Cohort) Genotype Preprocessing Pipeline
# =============================================================================
# Purpose: QC, HWE filtering, and format normalization for BBJ genotypes
# Input: PLINK binary format (bed/bim/fam)
# Output: QC-filtered PLINK binary with standardized VCF
# =============================================================================

set -euo pipefail

# Parse command line arguments
RAW_PREFIX="${1:?Error: RAW_PREFIX not provided}"
FASTA_REF="${2:?Error: FASTA_REF not provided}"
MIND="${3:-0.01}"
GENO="${4:-0.01}"
MAF="${5:-0.05}"
HWE="${6:-1e-6}"
THREADS="${7:-4}"
PLINK2_MEMORY_MB="${8:-60000}"

detected_mem_mb=""
if [[ -n "${SLURM_MEM_PER_NODE:-}" && "${SLURM_MEM_PER_NODE}" != "0" ]]; then
	detected_mem_mb="${SLURM_MEM_PER_NODE}"
elif [[ -n "${SLURM_MEM_PER_CPU:-}" && -n "${SLURM_CPUS_PER_TASK:-}" && "${SLURM_MEM_PER_CPU}" != "0" && "${SLURM_CPUS_PER_TASK}" != "0" ]]; then
	typeset -i _mem_per_cpu="${SLURM_MEM_PER_CPU}"
	typeset -i _cpus_per_task="${SLURM_CPUS_PER_TASK}"
	detected_mem_mb=$(( _mem_per_cpu * _cpus_per_task ))
fi

if [[ -n "${detected_mem_mb}" ]]; then
	# Keep a 2GB headroom to reduce cgroup OOM-kill risk
	typeset -i _mem_headroom_mb=2048
	typeset -i _safe_mem_mb=$(( detected_mem_mb - _mem_headroom_mb ))
	if (( _safe_mem_mb < 1024 )); then
		_safe_mem_mb=1024
	fi
	PLINK2_MEMORY_MB=${_safe_mem_mb}
fi

echo "[$(date)] Runtime resources: THREADS=${THREADS}, PLINK2_MEMORY_MB=${PLINK2_MEMORY_MB}"

# Derived variables
RAW_BED="${RAW_PREFIX}.bed"
RAW_BIM="${RAW_PREFIX}.bim"
RAW_FAM="${RAW_PREFIX}.fam"

# Check inputs
[ -f "$RAW_BED" ] || { echo "Error: $RAW_BED not found"; exit 1; }
[ -f "$RAW_BIM" ] || { echo "Error: $RAW_BIM not found"; exit 1; }
[ -f "$RAW_FAM" ] || { echo "Error: $RAW_FAM not found"; exit 1; }
[ -f "$FASTA_REF" ] || { echo "Error: $FASTA_REF not found"; exit 1; }

echo "[$(date)] Starting BBJ genotype preprocessing..."

# =============================================================================
# PHASE 1: INITIAL QUALITY CONTROL
# =============================================================================

echo "[$(date)] PHASE 1: Initial QC (mind, geno, autosome, maf)..."

# Step 1.1: Initial variant and sample QC
plink2 \
	--bfile "$RAW_PREFIX" \
	--mind "$MIND" \
	--geno "$GENO" \
	--autosome \
	--maf "$MAF" \
	--memory "$PLINK2_MEMORY_MB" \
	--make-bed \
	--out bbj.b38.auto.prep.filt \
	--threads "$THREADS"

# =============================================================================
# PHASE 2: HWE QUALITY CONTROL
# =============================================================================

echo "[$(date)] PHASE 2: HWE filtering..."

plink2 \
	--bfile bbj.b38.auto.prep.filt \
	--hwe "$HWE" \
	--memory "$PLINK2_MEMORY_MB" \
	--make-bed \
	--out bbj.b38.auto.prep.hwe \
	--threads "$THREADS"

# =============================================================================
# PHASE 3: FORMAT CONVERSION AND NORMALIZATION
# =============================================================================

echo "[$(date)] PHASE 3: Format conversion and VCF normalization..."

# Step 4.1: Export to VCF format
plink2 \
	--bfile bbj.b38.auto.prep.hwe \
	--export vcf bgz id-paste=iid \
	--memory "$PLINK2_MEMORY_MB" \
	--out bbj.b38.auto.prep.hwe \
	--threads "$THREADS"

# Step 4.2: Streaming VCF normalization pipeline
for i in {1..22}; do echo -e "${i}\tchr${i}"; done > chr_rename.txt
bcftools annotate --rename-chrs chr_rename.txt bbj.b38.auto.prep.hwe.vcf.gz -Ou | \
bcftools norm \
	--multiallelics -any \
	--fasta-ref "$FASTA_REF" \
	--check-ref s \
	--threads "$THREADS" \
	-Ou | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o bbj.b38.auto.prep.setid.vcf.gz
bcftools index --threads "$THREADS" -t bbj.b38.auto.prep.setid.vcf.gz

# Step 4.3: Convert normalized VCF back to PLINK binary format
plink2 \
	--vcf bbj.b38.auto.prep.setid.vcf.gz \
	--double-id \
	--memory "$PLINK2_MEMORY_MB" \
	--make-bed \
	--out bbj.b38.auto.prep.tmp \
	--threads "$THREADS"

# Step 4.4: Prefix sample IDs with "bbj_" for cohort distinction
awk '{print $1, $2, "bbj_"$2, "bbj_"$2}' bbj.b38.auto.prep.tmp.fam > bbj.b38.auto.prep.update_ids.txt

plink2 \
	--bfile bbj.b38.auto.prep.tmp \
	--update-ids bbj.b38.auto.prep.update_ids.txt \
	--memory "$PLINK2_MEMORY_MB" \
	--make-bed \
	--out bbj.b38.auto.prep \
	--threads "$THREADS"

# =============================================================================
# CLEANUP: Remove all intermediate files
# =============================================================================

echo "[$(date)] Cleaning up intermediate files..."

rm -f chr_rename.txt bbj.b38.auto.prep.update_ids.txt
rm -f bbj.b38.auto.prep.filt.bed bbj.b38.auto.prep.filt.bim bbj.b38.auto.prep.filt.fam bbj.b38.auto.prep.filt.log bbj.b38.auto.prep.filt.nosex bbj.b38.auto.prep.filt.smiss
rm -f bbj.b38.auto.prep.hwe.bed bbj.b38.auto.prep.hwe.bim bbj.b38.auto.prep.hwe.fam bbj.b38.auto.prep.hwe.log bbj.b38.auto.prep.hwe.nosex
rm -f bbj.b38.auto.prep.hwe.vcf.gz bbj.b38.auto.prep.hwe.vcf.gz.tbi
rm -f bbj.b38.auto.prep.tmp.bed bbj.b38.auto.prep.tmp.bim bbj.b38.auto.prep.tmp.fam bbj.b38.auto.prep.tmp.log bbj.b38.auto.prep.tmp.nosex

echo "[$(date)] BBJ genotype preprocessing completed!"
echo "[$(date)] Output files:"
echo "  - bbj.b38.auto.prep.bed"
echo "  - bbj.b38.auto.prep.bim"
echo "  - bbj.b38.auto.prep.fam"
echo "  - bbj.b38.auto.prep.setid.vcf.gz"
echo "  - bbj.b38.auto.prep.setid.vcf.gz.tbi"
