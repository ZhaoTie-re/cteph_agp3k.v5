
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
params.high_ld_regions       = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/info/high-LD-regions-hg38-GRCh38_modified.txt'

// -----------------------------------------------------------------------------
// BBJ Configuration
// -----------------------------------------------------------------------------
params.bbj_mind              = 0.01
params.bbj_geno              = 0.01
params.bbj_maf               = 0.05
params.bbj_hwe               = 1e-6
params.bbj_threads           = 16
params.bbj_reuse_outputs     = true

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
	def bbj_publish_dir = "${params.out_dir}/09_bbj_preprocess/bbj_raw_qc_norm_plink"
	def bbj_reuse_outputs = params.bbj_reuse_outputs
	"""
	export PATH=/home/b/b37974/:\$PATH
	source activate ${params.conda_env_activate}

	expected_sig="mind=${params.bbj_mind}|geno=${params.bbj_geno}|maf=${params.bbj_maf}|hwe=${params.bbj_hwe}"
	echo "\${expected_sig}" > bbj.preprocess.signature.expected.txt

	if [ "${bbj_reuse_outputs}" = "true" ] && \
	   [ -s "${bbj_publish_dir}/bbj.b38.auto.prep.bed" ] && \
	   [ -s "${bbj_publish_dir}/bbj.b38.auto.prep.bim" ] && \
	   [ -s "${bbj_publish_dir}/bbj.b38.auto.prep.fam" ] && \
	   [ -s "${bbj_publish_dir}/bbj.b38.auto.prep.setid.vcf.gz" ] && \
	   [ -s "${bbj_publish_dir}/bbj.b38.auto.prep.setid.vcf.gz.tbi" ] && \
	   [ -s "${bbj_publish_dir}/bbj.preprocess.signature.txt" ] && \
	   cmp -s bbj.preprocess.signature.expected.txt "${bbj_publish_dir}/bbj.preprocess.signature.txt"; then
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
	PREPARE_BBJ_GENOTYPE(ch_bbj_raw)

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
	RUN_PIHAT_QC(ch_sample_qc_plink, ch_fprune_in, ch_sample_qc_metrics)
}

