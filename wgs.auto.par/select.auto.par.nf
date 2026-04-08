
nextflow.enable.dsl = 2

// -----------------------------------------------------------------------------
// Parameters Configuration
// -----------------------------------------------------------------------------
params.vcf_dir                = '/LARGE1/gr10478/platform/JHRPv5/workspace/pipeline/output/annovar.v5/annovar_concat'
params.out_dir                = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/wgs.auto.par/results'
params.sif_dir                = '/home/b/b37974/simg'
params.script_dir             = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/wgs.auto.par/scripts'
params.sample_list            = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/info/cteph_agp3k.v5.ls'
params.sample_info            = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/info/ph_agp3k_combined_final_20251112.xlsx'
params.nagasaki_pipeline_path = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline'
params.fasta                  = "${params.nagasaki_pipeline_path}/data/hs38DH.fa"
params.gatk_sif               = "${params.sif_dir}/gatk_latest.sif"
params.bbj_raw_prefix         = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/BBJ_genome_b38/00.raw_data/NewOE13_Auto.id.b38'

// -----------------------------------------------------------------------------
// Sample Info Configuration
// -----------------------------------------------------------------------------
params.conda_env_activate     = 'cteph_geno_pro'
params.sample_id_col          = 'ID JHRPv5'
params.phenotype_col          = 'Outcome'
params.phenotype_case_value   = 'PH'
params.phenotype_ctrl_value   = 'AGP3K'
params.sex_col                = 'Sex'
params.sex_female_value       = 'F'
params.sex_male_value         = 'M'
params.target_dp_col		  = 'Target DP'
params.dp_col 			      = 'DP'

// -----------------------------------------------------------------------------
// QC Configuration
// -----------------------------------------------------------------------------
params.sample_qc_config       = "${params.script_dir}/sample_qc_config.json"
params.variant_qc_vmiss_config = "${params.script_dir}/vqc_config_vmiss.json"
params.variant_qc_hwe_config   = "${params.script_dir}/vqc_config_hwe.json"
params.high_ld_regions       = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/info/high-LD-regions-hg38-GRCh38_modified.txt'
// true: in RUN_VARIANT_QC, exclude IIDs with SELECTED_FOR_REMOVAL=true from PI_HAT vertex-cover TSV for HWE calculations only (VMISS and AAF always use full samples)
params.variant_qc_exclude_pihat_for_hwe = true

// -----------------------------------------------------------------------------
// BBJ Configuration
// -----------------------------------------------------------------------------
params.bbj_mind              = 0.01
params.bbj_geno              = 0.01
params.bbj_maf               = 0.05
params.bbj_hwe               = 1e-6
params.bbj_threads           = 16
params.bbj_reuse_outputs     = true

//-----------------------------------------------------------------------------
// PopGMM Configuration for ancestry inference
//------------------------------------------------------------------------------
params.popgmm = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/PopGMM_output/cluster2_highconf_fid_iid_conf_ge_0p95.tsv'

// -----------------------------------------------------------------------------
// Processes
// -----------------------------------------------------------------------------

process PREPARE_VCF {
	executor 'slurm'
	queue 'gr10478b'
	time '36h'
	tag "${chr}"

	publishDir "${params.out_dir}/01_prepare_vcf", mode: 'symlink'

	input:
	val chr
	path sample_list

	output:
	tuple val(chr), path("${chr}.selected.pass.norm_split.setid.vcf.gz"), path("${chr}.selected.pass.norm_split.setid.vcf.gz.tbi")

	script:
	def input_vcf = chr == 'PAR' ? "${params.vcf_dir}/annovar.v5.PAR.hg38_multianno.vcf.gz" : "${params.vcf_dir}/annovar.v5.${chr}.hg38_multianno.vcf.gz"
	def output_vcf = "${chr}.selected.pass.norm_split.setid.vcf.gz"

	"""
	# Step 1: Select target samples
	# Step 2: Retain only PASS variants
	# Step 3: Filter variants with MAC >= 1 (remove monomorphic sites, both all-ref and all-alt)
	# Step 4: Normalize and split multiallelic variants
	# Step 5: Remove spanning deletion alleles (where ALT is *)
	# Step 6: Set variant IDs to CHROM:POS:REF:ALT format
	bcftools view ${input_vcf} --threads 4 -S ${sample_list} --force-samples -Ou | \
		bcftools view --threads 4 -f"PASS" -Ou | \
		bcftools view --threads 4 --min-ac 1:minor -Ou | \
		bcftools norm \
			--multiallelics -any \
			--fasta-ref ${params.fasta} \
			--check-ref s \
			--threads 4 \
			-Ou | \
		bcftools filter --threads 4 -e 'ALT="*"' -Ou | \
		bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o ${output_vcf}

	bcftools index --threads 4 -t ${output_vcf}
	"""
}

process FILTER_VQC {
	executor 'slurm'
	queue 'gr10478b'
	time '36h'
	tag "${chr}"

	publishDir "${params.out_dir}/02_filter_vqc", mode: 'symlink'

	input:
	tuple val(chr), path(vcf), path(vcf_tbi)

	output:
	tuple val(chr), path("${chr}.vqc.vcf.gz"), path("${chr}.vqc.vcf.gz.tbi")

	script:
	def filtered_vcf = "${chr}.vqc.vcf.gz"

	"""
	# Filter variants based on quality metrics:
	# VQSLOD > 10: Variant Quality Score Log-Odds (confidence in variant call)
	# MQ > 58.75: Mapping Quality (alignment quality of reads supporting the variant)
	bcftools view ${vcf} --threads 4 -i 'VQSLOD > 10 & MQ > 58.75' -Oz -o ${filtered_vcf}
	bcftools index --threads 4 -t ${filtered_vcf}
	"""
}

process ANNOTATE_AF_NORM_GT {
	executor 'slurm'
	queue 'gr10478b'
	time '36h'
	tag "${chr}"

	publishDir "${params.out_dir}/03_annotate_af_norm_gt", mode: 'symlink'

	input:
	tuple val(chr), path(vcf), path(vcf_tbi)

	output:
	tuple val(chr), path("${chr}.vqc.af.gtnorm.vcf.gz"), path("${chr}.vqc.af.gtnorm.vcf.gz.tbi")

	script:
	def tmp_vcf = "${chr}.vqc.af.tmp.vcf.gz"
	def output_vcf = "${chr}.vqc.af.gtnorm.vcf.gz"

	"""
	# Step 1: Add AlleleFraction (old version: AlleleBalance) annotation using GATK VariantAnnotator
	singularity exec \
		--bind /LARGE0:/LARGE0 \
		--bind /LARGE1:/LARGE1 \
		${params.gatk_sif} gatk --java-options "-Xmx8G -XX:ParallelGCThreads=4" VariantAnnotator \
		-R ${params.fasta} \
		-V ${vcf} \
		-O ${tmp_vcf} \
		-A AlleleFraction \
		--create-output-variant-index true

	# Step 2: Unphase and sort all genotypes with bcftools +setGT
	bcftools +setGT ${tmp_vcf} -Ou -- -t a -n u | bcftools view --threads 4 -Oz -o ${output_vcf}
	bcftools index --threads 4 -t ${output_vcf}

	# Clean up temp files
	rm -f ${tmp_vcf}*
	"""
}

process FILTER_GENOTYPE {
	executor 'slurm'
	queue 'gr10478b'
	time '36h'
	tag "${chr}"

	publishDir "${params.out_dir}/04_filter_genotype", mode: 'symlink'

	input:
	tuple val(chr), path(vcf), path(tbi)

	output:
	tuple val(chr), path("${chr}.gt_qc.norm.vcf.gz"), path("${chr}.gt_qc.norm.vcf.gz.tbi")

	script:
	def tmp_vcf = "${chr}.gt_qc.tmp.vcf.gz"
	def final_vcf = "${chr}.gt_qc.norm.vcf.gz"

	"""
	# Step 1: Normalization (sed nan -> NaN) and re-indexing
	bcftools view ${vcf} --threads 4 | sed 's/nan/NaN/g' | bgzip > ${tmp_vcf}
	bcftools index --threads 4 -t ${tmp_vcf}

	# Step 2: GATK VariantFiltration for autosomes and PAR with unified thresholds
	singularity exec \
		--bind /LARGE0:/LARGE0 \
		--bind /LARGE1:/LARGE1 \
		${params.gatk_sif} gatk --java-options "-Xmx8G -XX:ParallelGCThreads=4" VariantFiltration \
		-R ${params.fasta} \
		-V ${tmp_vcf} \
		-O ${final_vcf} \
		--genotype-filter-name "LowGQ" \
		--genotype-filter-expression "GQ < 20" \
		--genotype-filter-name "LowDP" \
		--genotype-filter-expression "DP < 8" \
		--genotype-filter-name "ABB_outlier" \
		--genotype-filter-expression "isHet == 1 && (AF < 0.2 || AF > 0.8)" \
		--genotype-filter-name "ABB_NaN" \
		--genotype-filter-expression "AF == 'NaN'" \
		--set-filtered-genotype-to-no-call true \
		--create-output-variant-index true

	# Step 3: Delete tmp_vcf
	rm -f ${tmp_vcf}*
	"""
}

process VCF_TO_PLINK {
	executor 'slurm'
	queue 'gr10478b'
	time '12h'
	tag "${chr}"

	publishDir "${params.out_dir}/05_vcf_to_plink", mode: 'symlink'

	input:
	tuple val(chr), path(vcf), path(tbi)

	output:
	path("${chr}.plink.bed"), emit: bed
	path("${chr}.plink.bim"), emit: bim
	path("${chr}.plink.fam"), emit: fam

	script:
	def mac_vcf = "${chr}.mac_filtered.vcf.gz"
	def out_prefix = "${chr}.plink"
	def update_tsv = "${chr}.plink.update.tsv"
	def update_script = "${params.script_dir}/build_plink_update_table.py"
	def split_par_opt = chr == 'PAR' ? '--split-par b38' : ''

	"""
	export PATH=/home/b/b37974/:\$PATH
	source activate ${params.conda_env_activate}

	# 1. Filter out MAC < 1 variants after GT filtering
	bcftools view ${vcf} --threads 4 --min-ac 1:minor -Oz -o ${mac_vcf}

	# 2. Build sex/pheno update table from sample_info.xlsx
	python ${update_script} \
		--xlsx ${params.sample_info} \
		--out ${update_tsv} \
		--sample-id-col "${params.sample_id_col}" \
		--pheno-col "${params.phenotype_col}" \
		--case-value "${params.phenotype_case_value}" \
		--ctrl-value "${params.phenotype_ctrl_value}" \
		--sex-col "${params.sex_col}" \
		--female-value "${params.sex_female_value}" \
		--male-value "${params.sex_male_value}"

	# 3. Convert to PLINK format, update sex and phenotype in one step
	plink2 \
		--vcf ${mac_vcf} \
		--double-id \
		${split_par_opt} \
		--update-sex ${update_tsv} \
		--pheno ${update_tsv} \
		--pheno-name PHENO \
		--make-bed \
		--out ${out_prefix} \
		--threads 4

	# Clean up temp files
	rm -f ${mac_vcf} ${update_tsv}
	"""
}

process MERGE_PLINK {
	executor 'slurm'
	queue 'gr10478b'
	time '24h'

	publishDir "${params.out_dir}/06_merged_plink", mode: 'symlink'

	input:
	path beds
	path bims
	path fams

	output:
	tuple path("cteph_agp3k_v5_wgs_merged.bed"), path("cteph_agp3k_v5_wgs_merged.bim"), path("cteph_agp3k_v5_wgs_merged.fam")

	script:
	"""
	export PATH=/home/b/b37974/:\$PATH

	# Create merge list in chromosome order (chr1-22 only)
	rm -f merge_list.txt
	for i in {1..22}; do
		if [ -f "chr\${i}.plink.bed" ]; then
			echo "chr\${i}.plink" >> merge_list.txt
		fi
	done

	# Merge using plink2
	plink2 \
		--pmerge-list merge_list.txt bfile \
		--make-bed \
		--out cteph_agp3k_v5_wgs_merged \
		--threads 4
	"""
}

process BUILD_SAMPLE_QC_TABLE {
	executor 'slurm'
	queue 'gr10478b'
	time '12h'

	publishDir "${params.out_dir}/07_sample_qc/metrics", mode: 'symlink'

	input:
	tuple path(merged_bed), path(merged_bim), path(merged_fam)

	output:
	path("*.sample_qc_metrics.tsv")
	path("cteph_agp3k_v5_wgs_merged.Fprune.prune.in")

	script:
	def merged_prefix = "cteph_agp3k_v5_wgs_merged"
	def smiss_file = "${merged_prefix}.smiss"
	def het_file = "${merged_prefix}.het"
	def qc_table = "${merged_prefix}.sample_qc_metrics.tsv"
	def qc_script = "${params.script_dir}/build_sample_qc_metrics_table.py"

	"""
	export PATH=/home/b/b37974/:\$PATH
	source activate ${params.conda_env_activate}

	# 1. Compute sample missingness from merged PLINK files
	plink2 \
		--bfile ${merged_prefix} \
		--autosome \
		--missing sample-only \
		--out ${merged_prefix} \
		--threads 8

	# 2. Select variants for F calculation: remove high-LD regions and perform LD pruning (50 5 0.2)
	plink2 \
		--bfile ${merged_prefix} \
		--autosome \
		--snps-only just-acgt \
		--maf 0.05 \
		--exclude range ${params.high_ld_regions} \
		--indep-pairwise 50 5 0.2 \
		--out ${merged_prefix}.Fprune \
		--threads 8

	# 3. Compute per-sample heterozygosity coefficient (F) on pruned autosomal SNPs
	plink2 \
		--bfile ${merged_prefix} \
		--extract ${merged_prefix}.Fprune.prune.in \
		--autosome \
		--het \
		--out ${merged_prefix} \
		--threads 8

	# 4. Build sample-level QC metrics table
	python ${qc_script} \
		--xlsx ${params.sample_info} \
		--fam ${merged_prefix}.fam \
		--smiss ${smiss_file} \
		--het ${het_file} \
		--out ${qc_table} \
		--sample-id-col "${params.sample_id_col}" \
		--target-dp-col "${params.target_dp_col}" \
		--dp-col "${params.dp_col}"
	"""
}

process RUN_SAMPLE_QC_FROM_METRICS {
	executor 'slurm'
	queue 'gr10478b'
	time '12h'

	publishDir "${params.out_dir}/07_sample_qc/run_qc", mode: 'symlink'

	input:
	tuple path(merged_bed), path(merged_bim), path(merged_fam)
	path sample_qc_metrics_tsv

	output:
	path("cteph_agp3k_v5_wgs_merged.sample_qc.detail.tsv")
	path("cteph_agp3k_v5_wgs_merged.sample_qc.remove.id")
	path("cteph_agp3k_v5_wgs_merged.sample_qc.keep.id")
	path("cteph_agp3k_v5_wgs_merged.sample_qc.summary.json")
	path("cteph_agp3k_v5_wgs_merged.sample_qc.summary.txt")
	path("cteph_agp3k_v5_wgs_merged.sample_qc.png")
	tuple path("cteph_agp3k_v5_wgs_merged.sample_qc.bed"), path("cteph_agp3k_v5_wgs_merged.sample_qc.bim"), path("cteph_agp3k_v5_wgs_merged.sample_qc.fam")

	script:
	def merged_prefix = "cteph_agp3k_v5_wgs_merged"
	def qc_script = "${params.script_dir}/run_sample_qc.py"
	def out_prefix = "${merged_prefix}.sample_qc"
	def tmp_prefix = "${out_prefix}.tmp"

	"""
	export PATH=/home/b/b37974/:\$PATH
	source activate ${params.conda_env_activate}

	# 1. Apply sample QC rules to build remove/keep list and generate visualization
	python ${qc_script} \
		--metrics-tsv ${sample_qc_metrics_tsv} \
		--config-json ${params.sample_qc_config} \
		--out-prefix ${merged_prefix} \
		--sample-info-xlsx ${params.sample_info} \
		--sample-id-col "${params.sample_id_col}" \
		--phenotype-col "${params.phenotype_col}" \
		--case-value "${params.phenotype_case_value}" \
		--ctrl-value "${params.phenotype_ctrl_value}" \
		--case-label CTEPH

	# 2. Apply remove list to merged PLINK (sample-level QC)
	plink2 \
		--bfile ${merged_prefix} \
		--remove ${merged_prefix}.sample_qc.remove.id \
		--make-bed \
		--out ${tmp_prefix} \
		--threads 8

	# 3. Remove monomorphic variants after sample removal
	plink2 \
		--bfile ${tmp_prefix} \
		--mac 1 \
		--make-bed \
		--out ${out_prefix} \
		--threads 8

	# 4. Clean up temporary PLINK files
	rm -f ${tmp_prefix}.bed ${tmp_prefix}.bim ${tmp_prefix}.fam ${tmp_prefix}.log ${tmp_prefix}.nosex
	"""
}

process RUN_PIHAT_QC {
	executor 'slurm'
	queue 'gr10478b'
	time '12h'

	publishDir "${params.out_dir}/08_pi_hat_qc", mode: 'symlink'

	input:
	// Genotypes after sample-level QC from RUN_SAMPLE_QC_FROM_METRICS
	tuple path(qc_bed), path(qc_bim), path(qc_fam)
	// Pruned SNP list (Fprune.prune.in) used for PI_HAT calculation
	path prune_in
	// Sample-level QC metrics table providing SMISS etc.
	path sample_qc_metrics_tsv

	output:
	// High-PI_HAT pairs table
	path "*.pi_hat.pairs.tsv"
	// Per-sample vertex-cover annotation table
	path "*.pi_hat.vertex_cover_samples.tsv"
	// Human-readable log describing the strategy
	path "*.pi_hat.log.txt"
	// Network visualization (PNG)
	path "*.pi_hat.network.png"

	script:
	def bed_prefix = qc_bed.baseName
	def genome_prefix = "${bed_prefix}.pi_hat_genome"
	def qc_script = "${params.script_dir}/run_pihat_network_qc.py"

	"""
	export PATH=/home/b/b37974/:\$PATH
	source activate ${params.conda_env_activate}

	# 1. Compute pairwise PI_HAT using plink (not plink2),
	#    restricted to LD-pruned autosomal SNPs from BUILD_SAMPLE_QC_TABLE (Fprune.prune.in)
	plink \
		--bfile ${bed_prefix} \
		--extract ${prune_in} \
		--genome \
		--out ${genome_prefix} \
		--threads 8

	# 2. Build kinship network and weighted vertex cover summary
	python ${qc_script} \
		--genome ${genome_prefix}.genome \
		--metrics-tsv ${sample_qc_metrics_tsv} \
		--sample-info-xlsx ${params.sample_info} \
		--sample-id-col "${params.sample_id_col}" \
		--phenotype-col "${params.phenotype_col}" \
		--case-value "${params.phenotype_case_value}" \
		--ctrl-value "${params.phenotype_ctrl_value}" \
		--pi-hat-threshold 0.20 \
		--out-prefix ${bed_prefix}.pi_hat
	"""
}

process PREPARE_BBJ_GENOTYPE {
	executor 'slurm'
	queue 'gr10478b'
	time '24h'

	publishDir "${params.out_dir}/09_bbj_preprocess/bbj_raw_qc_norm_plink", mode: 'symlink'

	input:
	tuple path(raw_bed), path(raw_bim), path(raw_fam)

	output:
	tuple path("bbj.b38.auto.prep.bed"), path("bbj.b38.auto.prep.bim"), path("bbj.b38.auto.prep.fam")
	path("bbj.b38.auto.prep.setid.vcf.gz")
	path("bbj.b38.auto.prep.setid.vcf.gz.tbi")
	path("bbj.preprocess.signature.txt")

	script:
	def raw_prefix = raw_bed.baseName
	def bbj_script = "${params.script_dir}/prepare_bbj_genotype.sh"
	def signature_cmp_script = "${params.script_dir}/compare_bbj_preprocess_signature.py"
	def bbj_publish_dir = "${params.out_dir}/09_bbj_preprocess/bbj_raw_qc_norm_plink"
	def bbj_reuse_outputs = params.bbj_reuse_outputs
	"""
	export PATH=/home/b/b37974/:\$PATH
	source activate ${params.conda_env_activate}

	expected_sig="mind=${params.bbj_mind}|geno=${params.bbj_geno}|maf=${params.bbj_maf}|hwe=${params.bbj_hwe}"
	echo "\${expected_sig}" > bbj.preprocess.signature.expected.txt

	sig_match=false
	if [ "${bbj_reuse_outputs}" = "true" ] && \
	   [ -s "${bbj_publish_dir}/bbj.b38.auto.prep.bed" ] && \
	   [ -s "${bbj_publish_dir}/bbj.b38.auto.prep.bim" ] && \
	   [ -s "${bbj_publish_dir}/bbj.b38.auto.prep.fam" ] && \
	   [ -s "${bbj_publish_dir}/bbj.b38.auto.prep.setid.vcf.gz" ] && \
	   [ -s "${bbj_publish_dir}/bbj.b38.auto.prep.setid.vcf.gz.tbi" ] && \
	   [ -s "${bbj_publish_dir}/bbj.preprocess.signature.txt" ]; then
		if python ${signature_cmp_script} bbj.preprocess.signature.expected.txt "${bbj_publish_dir}/bbj.preprocess.signature.txt"
		then
			sig_match=true
		fi
	fi

	if [ "\${sig_match}" = "true" ]; then
		echo "[\$(date)] Found existing BBJ outputs with matching signature in ${bbj_publish_dir}, skip recomputation and reuse files."
		ln -sf "${bbj_publish_dir}/bbj.b38.auto.prep.bed" bbj.b38.auto.prep.bed
		ln -sf "${bbj_publish_dir}/bbj.b38.auto.prep.bim" bbj.b38.auto.prep.bim
		ln -sf "${bbj_publish_dir}/bbj.b38.auto.prep.fam" bbj.b38.auto.prep.fam
		ln -sf "${bbj_publish_dir}/bbj.b38.auto.prep.setid.vcf.gz" bbj.b38.auto.prep.setid.vcf.gz
		ln -sf "${bbj_publish_dir}/bbj.b38.auto.prep.setid.vcf.gz.tbi" bbj.b38.auto.prep.setid.vcf.gz.tbi
		ln -sf "${bbj_publish_dir}/bbj.preprocess.signature.txt" bbj.preprocess.signature.txt
	else
		echo "[\$(date)] Running BBJ preprocessing script..."

	# Run BBJ preprocessing pipeline script
	zsh ${bbj_script} \
		${raw_prefix} \
		${params.fasta} \
		${params.bbj_mind} \
		${params.bbj_geno} \
		${params.bbj_maf} \
		${params.bbj_hwe} \
		${params.bbj_threads}
		cp bbj.preprocess.signature.expected.txt bbj.preprocess.signature.txt
	fi
	rm -f bbj.preprocess.signature.expected.txt
	"""
}

process RUN_VARIANT_QC {
	executor 'slurm'
	queue 'gr10478b'
	time '24h'

	publishDir "${params.out_dir}/10_variant_qc", mode: 'symlink'

	input:
	tuple path(qc_bed), path(qc_bim), path(qc_fam)
	path pihat_vertex_cover_tsv

	output:
	path("*.variant_qc_summary.tsv")
	path("*.vmiss_pass_variants.tsv")
	path("*.hwe_pass_variants.tsv")
	path("*.pass_variants.tsv")
	tuple path("*.variant_qc.bed"), path("*.variant_qc.bim"), path("*.variant_qc.fam")
	path("*.vmiss.*.png")
	path("*.hwe.png")

	script:
	def bed_prefix = qc_bed.baseName
	def output_prefix = "${bed_prefix}.variant_qc"
	def variant_qc_script = "${params.script_dir}/variant_qc_pipeline.py"
	def pihat_exclude_arg = params.variant_qc_exclude_pihat_for_hwe ? "--pihat-vertex-cover-tsv-for-hwe ${pihat_vertex_cover_tsv}" : ""
	"""
	export PATH=/home/b/b37974/:\$PATH
	source activate ${params.conda_env_activate}

	python ${variant_qc_script} \
		--bed-prefix ${bed_prefix} \
		--out-prefix ${output_prefix} \
		--sample-info-xlsx ${params.sample_info} \
		--sample-id-col "${params.sample_id_col}" \
		--target-dp-col "${params.target_dp_col}" \
		--phenotype-col "${params.phenotype_col}" \
		--case-value "${params.phenotype_case_value}" \
		--ctrl-value "${params.phenotype_ctrl_value}" \
		--vmiss-config ${params.variant_qc_vmiss_config} \
		--vmiss-mode dp \
		--maf-group ctrl \
		--hwe-config ${params.variant_qc_hwe_config} \
		--script-path ${params.script_dir} \
		--tmpdir "${output_prefix}_tmp" \
		--threads 16 \
		${pihat_exclude_arg} \
		--no-stratify-by-maf
	"""
}


process PREPARE_BBJ_PCA_BASE {
	executor 'slurm'
	queue 'gr10478b'
	time '24h'

	publishDir "${params.out_dir}/09_bbj_preprocess/bbj_pca_base", mode: 'symlink'

	input:
	// BBJ reference genotype after preprocessing
	tuple path(bbj_bed), path(bbj_bim), path(bbj_fam)
	// Case-control genotype after variant QC
	tuple path(vqc_bed), path(vqc_bim), path(vqc_fam)

	output:
	// Keep outputs as separate channels to maximize resume compatibility.
	path("bbj.pca.intersect.snps")
	path("bbj.pca_base.eigenvec")
	path("bbj.pca_base.eigenval")
	path("bbj.pca_base.eigenvec.allele")
	path("bbj.pca_base.acount")

	script:
	def bbj_prefix = bbj_bed.baseName
	def vqc_prefix = vqc_bed.baseName
	"""
	export PATH=/home/b/b37974/:\$PATH
	source activate ${params.conda_env_activate}

	# 1. On BBJ genotype, remove high-LD regions and perform LD pruning (50 5 0.2)
	plink2 \
		--bfile ${bbj_prefix} \
		--autosome \
		--snps-only just-acgt \
		--exclude range ${params.high_ld_regions} \
		--indep-pairwise 50 5 0.2 \
		--out bbj.pca_prune \
		--threads 16

	# 2. From case-control variant QC genotype, export the list of variant IDs
	plink2 \
		--bfile ${vqc_prefix} \
		--write-snplist \
		--out vqc_all_snps \
		--threads 16

	# 3. Take the intersection of pruned BBJ SNPs and variant-QC SNPs
	sort -u bbj.pca_prune.prune.in > bbj.prune.sorted
	sort -u vqc_all_snps.snplist > vqc.snps.sorted
	comm -12 bbj.prune.sorted vqc.snps.sorted > bbj.pca.intersect.snps

	# 4. Run PCA on BBJ genotype using the intersected SNP set
	plink2 \
		--bfile ${bbj_prefix} \
		--extract bbj.pca.intersect.snps \
		--freq counts \
		--pca 20 allele-wts approx \
		--out bbj.pca_base \
		--threads 16
	"""
}


process PROJECT_ONTO_BBJ_PCS {
	executor 'slurm'
	queue 'gr10478b'
	time '24h'

	publishDir "${params.out_dir}/11_bbj_projection", mode: 'symlink'

	input:
	// PCA base outputs from PREPARE_BBJ_PCA_BASE (separate channels)
	path bbj_pca_snps
	path bbj_pca_evec
	path bbj_pca_eval
	path bbj_pca_evec_allele
	path bbj_pca_acount
	// BBJ reference genotype after preprocessing
	tuple path(bbj_bed), path(bbj_bim), path(bbj_fam)
	// Case-control genotype after variant QC
	tuple path(vqc_bed), path(vqc_bim), path(vqc_fam)

	output:
	// Projected PCs (scores) for merged BBJ + case/control genotypes
	path("*.bbjproj.sscore")
	path("*.bbjproj.sscore.vars")
	// Figure 1: explained ratio / cumulative explained ratio
	path("*.bbjproj.variance_summary.png")
	// Figure 2: pairwise PC scatter plots
	path("*.bbjproj.pc_pairs.pdf")

	script:
	def bbj_prefix  = bbj_bed.baseName
	def vqc_prefix  = vqc_bed.baseName
	def proj_prefix = "${vqc_prefix}.bbjproj"
	def projection_plot_script = "${params.script_dir}/plot_bbj_projection.py"
	"""
	export PATH=/home/b/b37974/:\$PATH
	source activate ${params.conda_env_activate}

	# 1) Restrict BBJ and case-control genotypes to the common SNP set used by BBJ PCA base.
	plink2 \
		--bfile ${bbj_prefix} \
		--extract ${bbj_pca_snps} \
		--make-bed \
		--out bbj.proj.base \
		--threads 16

	plink2 \
		--bfile ${vqc_prefix} \
		--extract ${bbj_pca_snps} \
		--make-bed \
		--out vqc.proj.base \
		--threads 16

	# 2) Merge BBJ and case-control genotypes with plink1.9 --bmerge.
	#    This path is more stable than plink2 --pmerge-list for this use-case.
	plink \
		--bfile bbj.proj.base \
		--bmerge vqc.proj.base.bed vqc.proj.base.bim vqc.proj.base.fam \
		--make-bed \
		--keep-allele-order \
		--out merged.proj.base \
		--threads 16

	# 3) True projection onto BBJ PCs using allele weights and reference frequencies.
	plink2 \
		--bfile merged.proj.base \
		--read-freq ${bbj_pca_acount} \
		--score ${bbj_pca_evec_allele} 2 6 header-read no-mean-imputation variance-standardize list-variants \
		--score-col-nums 7-26 \
		--out ${proj_prefix} \
		--threads 16

	# Publication-style figures for BBJ PCA and projection results.
	python ${projection_plot_script} \
		--bbj-eigenval ${bbj_pca_eval} \
		--projected-sscore ${proj_prefix}.sscore \
		--sample-info ${params.sample_info} \
		--sample-id-col "${params.sample_id_col}" \
		--phenotype-col "${params.phenotype_col}" \
		--phenotype-case-value "${params.phenotype_case_value}" \
		--phenotype-ctrl-value "${params.phenotype_ctrl_value}" \
		--bbj-id-prefix "bbj_" \
		--bbj-label "BBJ" \
		--case-label "CTEPH" \
		--ctrl-label "AGP3K" \
		--max-pcs 20 \
		--out-prefix ${proj_prefix}
	"""
}


process POPGMM_SUBSET_AND_PLOT_BBJ_PROJECTION {
	executor 'slurm'
	queue 'gr10478b'
	time '24h'

	publishDir "${params.out_dir}/12_popgmm_subset_projection", mode: 'symlink'

	input:
	// Variant-QC genotype (source for PopGMM subsetting)
	tuple path(vqc_bed), path(vqc_bim), path(vqc_fam)
	// PopGMM high-confidence keep list (FID IID, no header)
	path popgmm_keep
	// Existing projection score file from PROJECT_ONTO_BBJ_PCS
	path projected_sscore
	// BBJ eigenvalues for variance/explained-ratio plotting
	path bbj_pca_eval

	output:
	path("*.popgmm.bed")
	path("*.popgmm.bim")
	path("*.popgmm.fam")
	path("*.popgmm.subset.log.txt")
	path("*.bbjproj.popgmm.variance_summary.png")
	path("*.bbjproj.popgmm.pc_pairs.pdf")

	script:
	def vqc_prefix = vqc_bed.baseName
	def subset_pre = "${vqc_prefix}.popgmm.keep"
	def subset_out = "${vqc_prefix}.popgmm"
	def proj_prefix = "${vqc_prefix}.bbjproj.popgmm"
	def subset_log = "${vqc_prefix}.popgmm.subset.log.txt"
	def projection_plot_script = "${params.script_dir}/plot_bbj_projection.py"
	"""
	export PATH=/home/b/b37974/:\$PATH
	source activate ${params.conda_env_activate}

	# 1) Subset variant-QC genotype by PopGMM samples
	plink2 \
		--bfile ${vqc_prefix} \
		--keep ${popgmm_keep} \
		--make-bed \
		--out ${subset_pre} \
		--threads 16

	# 2) Remove monomorphic variants after sample subset
	plink2 \
		--bfile ${subset_pre} \
		--mac 1 \
		--make-bed \
		--out ${subset_out} \
		--threads 16

	# 3) Log detailed sample/variant count changes
	before_samples=\$(wc -l < ${vqc_prefix}.fam)
	before_variants=\$(wc -l < ${vqc_prefix}.bim)
	after_keep_samples=\$(wc -l < ${subset_pre}.fam)
	after_keep_variants=\$(wc -l < ${subset_pre}.bim)
	after_mac_samples=\$(wc -l < ${subset_out}.fam)
	after_mac_variants=\$(wc -l < ${subset_out}.bim)

	{
		echo "[\$(date)] POPGMM subset + monomorphic-variant removal summary"
		echo "INPUT_BFILE_PREFIX: ${vqc_prefix}"
		echo "POPGMM_KEEP_FILE: ${popgmm_keep}"
		echo "STEP1_KEEP_PREFIX: ${subset_pre}"
		echo "STEP2_FINAL_PREFIX: ${subset_out}"
		echo ""
		echo "Counts (samples / variants):"
		echo "  Before PopGMM keep            : \${before_samples} / \${before_variants}"
		echo "  After PopGMM keep             : \${after_keep_samples} / \${after_keep_variants}"
		echo "  After monomorphic rm (--mac 1): \${after_mac_samples} / \${after_mac_variants}"
		echo ""
		echo "Delta (after - before):"
		echo "  Keep step   sample delta: \$((after_keep_samples - before_samples))"
		echo "  Keep step  variant delta: \$((after_keep_variants - before_variants))"
		echo "  MAC step    sample delta: \$((after_mac_samples - after_keep_samples))"
		echo "  MAC step   variant delta: \$((after_mac_variants - after_keep_variants))"
		echo "  Total       sample delta: \$((after_mac_samples - before_samples))"
		echo "  Total      variant delta: \$((after_mac_variants - before_variants))"
	} > ${subset_log}

	# 4) Replot BBJ projection with optional non-BBJ PopGMM filtering
	python ${projection_plot_script} \
		--bbj-eigenval ${bbj_pca_eval} \
		--projected-sscore ${projected_sscore} \
		--sample-info ${params.sample_info} \
		--sample-id-col "${params.sample_id_col}" \
		--phenotype-col "${params.phenotype_col}" \
		--phenotype-case-value "${params.phenotype_case_value}" \
		--phenotype-ctrl-value "${params.phenotype_ctrl_value}" \
		--bbj-id-prefix "bbj_" \
		--bbj-label "BBJ" \
		--case-label "CTEPH" \
		--ctrl-label "AGP3K" \
		--max-pcs 20 \
		--keep-non-bbj-iids ${popgmm_keep} \
		--out-prefix ${proj_prefix}

	# Clean temporary keep-only files
	rm -f ${subset_pre}.bed ${subset_pre}.bim ${subset_pre}.fam ${subset_pre}.log ${subset_pre}.nosex
	"""
}





// -----------------------------------------------------------------------------
// Workflow Execution
// -----------------------------------------------------------------------------

workflow {
	// 1. Initialize input channels: chr1-22 and PAR
	ch_chrs = channel.fromList((1..22).collect { chrNum -> "chr${chrNum}" } + ['PAR'])
	ch_sample_list = file(params.sample_list, checkIfExists: true)
	ch_bbj_raw = channel.of([
		file("${params.bbj_raw_prefix}.bed", checkIfExists: true),
		file("${params.bbj_raw_prefix}.bim", checkIfExists: true),
		file("${params.bbj_raw_prefix}.fam", checkIfExists: true)
	])

	// [BBJ-PREP] Raw genotype preprocessing for projection (outside main numbered flow)
	ch_bbj_prepped = PREPARE_BBJ_GENOTYPE(ch_bbj_raw)
	ch_bbj_plink   = ch_bbj_prepped[0]

	// 2. Execute process flow
	ch_prepared         = PREPARE_VCF(ch_chrs, ch_sample_list)
	ch_filtered_vqc     = FILTER_VQC(ch_prepared)
	ch_annotated_gtnorm = ANNOTATE_AF_NORM_GT(ch_filtered_vqc)
	ch_final_genotype   = FILTER_GENOTYPE(ch_annotated_gtnorm)

	// 3. Filter MAC and convert to PLINK format
	ch_plink_files      = VCF_TO_PLINK(ch_final_genotype)

	// 4. Merge all PLINK files (only chr1-22,now)
	ch_merged = MERGE_PLINK(
		ch_plink_files.bed.collect(),
		ch_plink_files.bim.collect(),
		ch_plink_files.fam.collect()
	)

	// 5. Build sample-level QC metrics table from merged PLINK outputs
	//    BUILD_SAMPLE_QC_TABLE has two outputs: metrics TSV and Fprune.prune.in
	ch_build_qc_all      = BUILD_SAMPLE_QC_TABLE(ch_merged)
	ch_sample_qc_metrics = ch_build_qc_all[0]
	ch_fprune_in         = ch_build_qc_all[1]

	// 6. Apply sample QC rules and generate post-QC PLINK files
	ch_sample_qc_all   = RUN_SAMPLE_QC_FROM_METRICS(ch_merged, ch_sample_qc_metrics)
	ch_sample_qc_plink = ch_sample_qc_all[6]

	// 7. Run PI_HAT-based relatedness QC on post-QC genotypes (annotation only)
	ch_pihat_all    = RUN_PIHAT_QC(ch_sample_qc_plink, ch_fprune_in, ch_sample_qc_metrics)
	ch_pihat_vertex = ch_pihat_all[1]

	// 8. Run variant QC on post-sample-QC genotypes
	ch_variant_qc_all   = RUN_VARIANT_QC(ch_sample_qc_plink, ch_pihat_vertex)
	ch_variant_qc_plink = ch_variant_qc_all[4]

	// 9. Build BBJ PCA base for future projection using intersected pruned SNPs
	ch_bbj_pca_base = PREPARE_BBJ_PCA_BASE(ch_bbj_plink, ch_variant_qc_plink)

	// 10. Project case/control samples onto PCs using the same SNP intersection
	ch_projected_all = PROJECT_ONTO_BBJ_PCS(
		ch_bbj_pca_base[0],
		ch_bbj_pca_base[1],
		ch_bbj_pca_base[2],
		ch_bbj_pca_base[3],
		ch_bbj_pca_base[4],
		ch_bbj_plink,
		ch_variant_qc_plink
	)

	// 11. PopGMM subset on variant-QC genotype + PopGMM-filtered replot from existing projection sscore
	ch_popgmm_keep = file(params.popgmm, checkIfExists: true)
	POPGMM_SUBSET_AND_PLOT_BBJ_PROJECTION(
		ch_variant_qc_plink,
		ch_popgmm_keep,
		ch_projected_all[0],
		ch_bbj_pca_base[2]
	)
}

