nextflow.enable.dsl = 2

// -----------------------------------------------------------------------------
// Runtime Configuration
// -----------------------------------------------------------------------------
params.genotypeRoot = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/wgs.auto.par/results'
params.tommoDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN'
params.cohortVcfSource = '/LARGE1/gr10478/platform/JHRPv4/workspace/pipeline/output/VQSR.v4'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/analysis/assoc_plink2/scripts'
params.resultsDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/analysis/assoc_plink2/results'
params.condaEnvActivate = 'cteph_geno_pro'

// -----------------------------------------------------------------------------
// Analysis Configuration
// -----------------------------------------------------------------------------
params.models = ['additive', 'dominant', 'recessive']
params.summaryPThreshold = 5e-8
params.plink2Threads = 16
params.phenoName = 'PHENO1'
params.covarNames = 'SEX,AGE_Z,PC1_AVG-PC7_AVG'
params.covarPcSource = 'bbj'   // own | bbj
params.onlySnpAssoc = false
params.includeAgeZ = params.covarNames
    .split(',')
    .collect { String s -> s.trim().toUpperCase() }
    .contains('AGE_Z')

// -----------------------------------------------------------------------------
// Relative Input Paths
// -----------------------------------------------------------------------------
params.bedPrefixRelPath = '14_fixed_model_prep/cteph_agp3k_v5_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.maf_ge_threshold'
params.phenoRelPath = '16_cov_pheno_prep/popgmm_subset_on_bbj_pcs.pheno.tsv'
params.covarRelPath = params.covarPcSource == 'bbj' \
    ? (params.includeAgeZ \
        ? '16_cov_pheno_prep/popgmm_subset_on_bbj_pcs.cov.sex_age_agez.tsv' \
        : '16_cov_pheno_prep/popgmm_subset_on_bbj_pcs.cov.sex.tsv') \
    : (params.includeAgeZ \
        ? '16_cov_pheno_prep/popgmm_relatedness_aware_projection.cov.sex_age_agez.tsv' \
        : '16_cov_pheno_prep/popgmm_relatedness_aware_projection.cov.sex.tsv')
params.ageNaRemoveRelPath = '16_cov_pheno_prep/popgmm_subset_on_bbj_pcs.age_na.fid_iid'
params.tommoVcfName = 'tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz'

// -----------------------------------------------------------------------------
// Process: Model Association Runs
// -----------------------------------------------------------------------------
process RunAssocPlink2Model {
    executor 'slurm'
    queue 'gr10478b'
    time '6d'
    tag "assoc_plink2: ${modelName}"

    publishDir "${params.resultsDir}/01.assoc_result/${modelName}", mode: 'symlink'

    input:
    tuple val(genotypeRoot), val(modelName)

    output:
    path('*.log'), emit: logs
    tuple val(modelName), path('*.glm.logistic*'), emit: glm

    script:
    def bedPrefix = "${genotypeRoot}/${params.bedPrefixRelPath}"
    def phenoFile = "${genotypeRoot}/${params.phenoRelPath}"
    def covarFile = "${genotypeRoot}/${params.covarRelPath}"
    def modelOpt = modelName == 'additive' ? '--glm' : "--glm ${modelName}"
    def snpOnlyOpt = params.onlySnpAssoc ? '--snps-only just-acgt' : ''
    def ageNaRemoveFile = "${genotypeRoot}/${params.ageNaRemoveRelPath}"
    def ageNaRemoveOpt = params.includeAgeZ ? "--remove ${ageNaRemoveFile}" : ''

    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.condaEnvActivate}
    covarTag=\$(python ${params.scriptDir}/covar_tag_from_names.py --covar-names "${params.covarNames}")
    outPrefix="cteph_agp3k.\${covarTag}.${modelName}"
    plink2 \\
        --bfile ${bedPrefix} \\
        --pheno ${phenoFile} \\
        --pheno-name ${params.phenoName} \\
        --covar ${covarFile} \\
        --covar-name "${params.covarNames}" \\
        ${ageNaRemoveOpt} \\
        ${snpOnlyOpt} \\
        ${modelOpt} omit-ref no-firth hide-covar \\
        --out \${outPrefix} \\
        --ci 0.95 \\
        --threads ${params.plink2Threads}
    """
}

// -----------------------------------------------------------------------------
// Process: Summary Integration
// -----------------------------------------------------------------------------
process BuildAssocPlink2Summary {
    executor 'slurm'
    queue 'gr10478b'
    time '6d'
    tag 'assoc_plink2: summary'

    publishDir "${params.resultsDir}/02.summary_vis", mode: 'symlink'

    input:
    tuple val(addModel), file(addResult)
    tuple val(domModel), file(domResult)
    tuple val(recModel), file(recResult)
    val(genotypeRoot)
    val(tommoDir)
    val(cohortVcfSource)

    output:
    path('*.csv')

    script:
    def bedPrefix = "${genotypeRoot}/${params.bedPrefixRelPath}"
    def tommoVcf = "${tommoDir}/${params.tommoVcfName}"

    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.condaEnvActivate}
    python ${params.scriptDir}/gwas_model_summary.py \\
        --bed-prefix ${bedPrefix} \\
        --tommo-vcf-file ${tommoVcf} \\
        --cohort-vcf-path ${cohortVcfSource} \\
        --p-threshold ${params.summaryPThreshold} \\
        --plink2-threads ${params.plink2Threads} \\
        --add-path ${addResult} \\
        --dom-path ${domResult} \\
        --rec-path ${recResult}
    """
}


// -----------------------------------------------------------------------------
// Process: PLINK2 Manhattan/QQ Visualization
// -----------------------------------------------------------------------------
process PlotAssocPlink2ModelResult {
    executor 'slurm'
    queue 'gr10478b'
    time '6d'
    tag "assoc_plink2: plot ${modelName}"

    publishDir "${params.resultsDir}/03.manhattan_qq/${modelName}", mode: 'symlink'

    input:
    tuple val(modelName), file(glmResult)

    output:
    tuple val(modelName), path('*.png'), emit: png
    tuple val(modelName), path('*.tsv'), emit: stats

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.condaEnvActivate}
    covarTag=\$(python ${params.scriptDir}/covar_tag_from_names.py --covar-names "${params.covarNames}")
    outPrefix="cteph_agp3k.\${covarTag}.${modelName}"
    python ${params.scriptDir}/plink2_manhattan_qq.py \
        --input ${glmResult} \
        --output-prefix \${outPrefix} \
        --title "PLINK2 ${modelName} model"
    """
}


// -----------------------------------------------------------------------------
// Workflow Execution
// -----------------------------------------------------------------------------
workflow {
    modelCh = channel.fromList(params.models)
    genotypeRootCh = channel.value(params.genotypeRoot)
    tommoDirCh = channel.value(params.tommoDir)
    cohortVcfSourceCh = channel.value(params.cohortVcfSource)

    modelInputCh = genotypeRootCh.combine(modelCh)

    RunAssocPlink2Model(modelInputCh)
    PlotAssocPlink2ModelResult(RunAssocPlink2Model.out.glm)

    addResultCh = RunAssocPlink2Model.out.glm.filter { row -> row[0] == 'additive' }
    domResultCh = RunAssocPlink2Model.out.glm.filter { row -> row[0] == 'dominant' }
    recResultCh = RunAssocPlink2Model.out.glm.filter { row -> row[0] == 'recessive' }

    BuildAssocPlink2Summary(
        addResultCh,
        domResultCh,
        recResultCh,
        genotypeRootCh,
        tommoDirCh,
        cohortVcfSourceCh
    )
}
