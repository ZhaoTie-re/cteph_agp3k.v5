nextflow.enable.dsl = 2

/*
 * SAIGE Association Pipeline (Nextflow DSL2)
 *
 * Overview:
 * 1) Prepare shared inputs (LD pruning)
 * 2) Full-GRM branch (mandatory): BGEN prep -> null model -> chr-wise association -> merge -> plots
 * 3) Sparse-GRM branch (optional): sparse GRM -> null model -> chr-wise association -> merge -> plots
 *
 * Naming convention:
 *   ${projectName}.${covarTag}.${grmType}.*
 * Example:
 *   cteph_agp3k.v5.saige.sex.10pc.fullGrm.chr1.assoc.txt
 *
 * Output layout:
 *   ${resultsDir}/00.prep
 *   ${resultsDir}/01.fullGrm
 *   ${resultsDir}/02.sparseGrm
 */

// -----------------------------------------------------------------------------
// Runtime Configuration
// -----------------------------------------------------------------------------
// Absolute paths and runtime environment.
params.genotypeRoot = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/wgs.auto.par/results'
params.infoDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/info'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/analysis/assoc_saige/scripts'
params.resultsDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/analysis/assoc_saige/results'
params.condaEnvActivate = 'saige'

// -----------------------------------------------------------------------------
// Analysis Configuration
// -----------------------------------------------------------------------------
// Core analysis controls.
params.projectName = 'cteph_agp3k.v5.saige'
params.phenoName = 'PHENO1'
params.covarNames = 'SEX,AGE_Z,PC1_AVG,PC2_AVG,PC3_AVG,PC4_AVG,PC5_AVG,PC6_AVG,PC7_AVG'
params.covarPcSource = 'bbj' // own | bbj
params.runSparse = true // true: run sparse GRM branch; false: full GRM only
params.includeAgeZ = params.covarNames
    .split(',')
    .collect { String s -> s.trim().toUpperCase() }
    .contains('AGE_Z')

// -----------------------------------------------------------------------------
// Input Path Configuration
// -----------------------------------------------------------------------------
// Relative input structure under genotypeRoot/infoDir.
params.randomModelDirRelPath = '15_random_model_prep'
params.covPhenoDirRelPath = '16_cov_pheno_prep'

params.plinkPrefixName = 'cteph_agp3k_v5_wgs_merged.sample_qc.variant_qc.popgmm.random_model'
params.phenoFileName = 'popgmm_subset_on_bbj_pcs.pheno.tsv'
params.covFileName = params.covarPcSource == 'bbj' \
    ? (params.includeAgeZ \
        ? 'popgmm_subset_on_bbj_pcs.cov.sex_age_agez.tsv' \
        : 'popgmm_subset_on_bbj_pcs.cov.sex.tsv') \
    : (params.includeAgeZ \
        ? 'popgmm_relatedness_aware_projection.cov.sex_age_agez.tsv' \
        : 'popgmm_relatedness_aware_projection.cov.sex.tsv')
params.ageNaRemoveFileName = 'popgmm_subset_on_bbj_pcs.age_na.fid_iid'
params.highLdFileName = 'high-LD-regions-hg38-GRCh38_modified.txt'


/*
 * Stage 0: Shared Preprocessing
 * Purpose: Build LD-pruned marker list used by both GRM branches.
 */
process LdPruning {
    executor 'slurm'
    queue 'gr10478b'
    time '6h'
    publishDir "${params.resultsDir}/00.prep/01.ldPruning", mode: 'symlink'

    input:
    tuple val(prefix), path(bed), path(bim), path(fam)
    path(high_ld)

    output:
    path('*.prune.in'), emit: pruneIn

    script:
    def outPrefix = "${prefix}.ldPruning"
    def normalizeHighLdScript = "${params.scriptDir}/normalize_high_ld_chr_prefix.sh"
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.condaEnvActivate}
    set -euo pipefail

    zsh ${normalizeHighLdScript} ${bim} ${high_ld} high_ld_regions.txt

    # Generate LD-pruned marker set
    
    plink2 \
        --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --exclude range high_ld_regions.txt \
        --maf 0.01 \
        --snps-only just-acgt \
        --indep-pairwise 50 10 0.2 \
        --out ${outPrefix} \
        --threads 8
    """
}

/*
 * Stage 0b (Shared, conditional): Remove AGE_Z-missing samples once
 * from the base PLINK dataset so all downstream steps use a consistent sample set.
 */
process RemoveAgeNaSamples {
    executor 'slurm'
    queue 'gr10478b'
    time '6h'
    publishDir "${params.resultsDir}/00.prep/00.ageNaRemove", mode: 'symlink'

    input:
    tuple val(prefix), path(bed), path(bim), path(fam)
    path(age_na_remove_file)

    output:
    tuple val(prefix), path('plink.agezrm.bed'), path('plink.agezrm.bim'), path('plink.agezrm.fam'), emit: filteredPlink

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.condaEnvActivate}
    set -euo pipefail

    plink2 \
        --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --remove ${age_na_remove_file} \
        --make-bed \
        --out plink.agezrm \
        --threads 8
    """
}

/*
 * Stage 1 (Full): Prepare per-chromosome BGEN + BGI + sample files.
 */
process SplitPlinkToBgen {
    executor 'slurm'
    queue 'gr10478b'
    time '6h'
    publishDir "${params.resultsDir}/01.fullGrm/01.preparedBgen", mode: 'symlink'

    input:
    tuple val(prefix), path(bed), path(bim), path(fam)
    each chr

    output:
    tuple val(chr), path('*.bgen'), path('*.bgen.bgi'), path('*.sample'), emit: bgenFiles

    script:
    def covarTagScript = "${params.scriptDir}/covar_tag_from_names.py"
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.condaEnvActivate}
    set -euo pipefail
    covarTag=\$(python3 ${covarTagScript} --covar-names "${params.covarNames}")
    outFile="${params.projectName}.\${covarTag}.fullGrm.chr${chr}"

    # Split by chromosome and export BGEN for SAIGE Step 2
    # bgen-1.2 bits=8: compatible and compact dosage encoding for SAIGE Step 2
    # ref-first: keep allele ordering consistent with SAIGE full-GRM settings
    # id-paste=iid: generate stable per-sample IDs aligned with IID-based phenotype/covariate tables
    
    plink2 \
        --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --chr ${chr} \
        --export bgen-1.2 bits=8 ref-first id-paste=iid \
        --out \$outFile \
        --threads 4
        
    # Create BGEN index for random access in SAIGE Step 2
    bgenix -g \$outFile.bgen -index -clobber
    """
}

/*
 * Stage 2 (Full): Fit SAIGE null model with full GRM.
 */
process FitNullGlmmFullGrm {
    executor 'slurm'
    queue 'gr10478b'
    time '72h'
    publishDir "${params.resultsDir}/01.fullGrm/02.nullModel", mode: 'symlink'

    input:
    tuple val(prefix), path(bed), path(bim), path(fam)
    path(prune_in)
    path(pheno_file)
    path(cov_file)
    path(age_na_remove_file)

    output:
    tuple path('*.rda'), path('*.varianceRatio.txt'), emit: nullModel

    script:
    def mergeScript = "${params.scriptDir}/merge_pheno_cov.py"
    def covarTagScript = "${params.scriptDir}/covar_tag_from_names.py"
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.condaEnvActivate}
    set -euo pipefail
    covarTag=\$(python3 ${covarTagScript} --covar-names "${params.covarNames}")
    fullNullPrefix="${params.projectName}.\${covarTag}.fullGrm.null"

    # Build temporary pruned PLINK dataset
    plink2 \
        --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --extract ${prune_in} \
        --make-bed \
        --out ${prefix}.pruned \
        --threads 8

    # Merge phenotype and covariates once for Step 1
    python3 ${mergeScript} \
        --pheno "${pheno_file}" \\
        --cov "${cov_file}" \\
        --pheno_col "${params.phenoName}" \\
        --cov_list "${params.covarNames}" \
        --sex_col "SEX" \\
        --out "merged_pheno_cov.txt"

    if [[ "${params.includeAgeZ}" == "true" ]]; then
        awk 'NR==FNR{rm[\$2]=1; next} FNR==1 || !(\$2 in rm)' ${age_na_remove_file} merged_pheno_cov.txt > merged_pheno_cov.filtered.txt
    else
        cp merged_pheno_cov.txt merged_pheno_cov.filtered.txt
    fi

    # SAIGE Step 1 (full GRM)
    
    step1_fitNULLGLMM.R \
        --plinkFile=${prefix}.pruned \
        --phenoFile=merged_pheno_cov.filtered.txt \
        --phenoCol=${params.phenoName} \
        --covarColList=${params.covarNames} \
        --sexCol=SEX \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --outputPrefix=\${fullNullPrefix} \
        --nThreads=32 \
        --isDiagofKinSetAsOne=FALSE \
        --numRandomMarkerforVarianceRatio=200 \
        --skipVarianceRatioEstimation=FALSE \
        --useSparseGRMtoFitNULL=FALSE \
        --IsOverwriteVarianceRatioFile=TRUE \
        --isCovariateOffset=FALSE
    """
}

/*
 * Stage 3 (Full): Run per-chromosome association tests on BGEN.
 */
process AssocTestFullGrm {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    publishDir "${params.resultsDir}/01.fullGrm/03.assocByChr", mode: 'symlink'

    input:
    tuple path(model_file), path(variance_ratio), val(chr), path(bgen_file), path(bgen_index), path(sample_file)

    output:
    tuple val(chr), path('*.assoc.txt'), emit: assocResults

    script:
    def covarTagScript = "${params.scriptDir}/covar_tag_from_names.py"
    
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.condaEnvActivate}
    set -euo pipefail
    covarTag=\$(python3 ${covarTagScript} --covar-names "${params.covarNames}")
    outFile="${params.projectName}.\${covarTag}.fullGrm.chr${chr}.assoc.txt"
    
    # Convert PLINK .sample to headerless sample ID list expected by SAIGE
    
    awk 'NR>2 {print \$1}' ${sample_file} > sample_ids.txt
    
    # SAIGE Step 2 (full GRM, BGEN input)
    
    step2_SPAtests.R \
        --bgenFile=${bgen_file} \
        --bgenFileIndex=${bgen_index} \
        --sampleFile=sample_ids.txt \
        --AlleleOrder=ref-first \
        --SAIGEOutputFile=\$outFile \
        --chrom=${chr} \
        --GMMATmodelFile=${model_file} \
        --varianceRatioFile=${variance_ratio} \
        --is_Firth_beta=FALSE \
        --LOCO=TRUE \
        --is_output_moreDetails=TRUE
    """
}

/*
 * Stage 1 (Sparse): Build sparse GRM matrix from pruned markers.
 */
process CalcSparseGrm {
    executor 'slurm'
    queue 'gr10478b'
    time '48h'
    publishDir "${params.resultsDir}/02.sparseGrm/01.sparseGrm", mode: 'symlink'

    input:
    tuple val(prefix), path(bed), path(bim), path(fam)
    path(prune_in)

    output:
    tuple path('*.sparseGRM.mtx'), path('*.sparseGRM.mtx.sampleIDs.txt'), emit: sparseGrm

    script:
    def outPrefix = prefix
    def covarTagScript = "${params.scriptDir}/covar_tag_from_names.py"
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.condaEnvActivate}
    set -euo pipefail
    covarTag=\$(python3 ${covarTagScript} --covar-names "${params.covarNames}")
    sparsePrefix="${params.projectName}.\${covarTag}.sparseGrm"

    # Build temporary pruned PLINK dataset
    plink2 \
        --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --extract ${prune_in} \
        --make-bed \
        --out ${outPrefix}.pruned \
        --threads 8

    # Compute sparse GRM
    createSparseGRM.R \
        --plinkFile=${outPrefix}.pruned \
        --nThreads=8 \
        --outputPrefix=\${sparsePrefix} \
        --numRandomMarkerforSparseKin=5000 \
        --relatednessCutoff=0.05 \
        --isDiagofKinSetAsOne=FALSE
    """
}

/*
 * Stage 2 (Sparse): Fit SAIGE null model with sparse GRM.
 */
process FitNullGlmmSparseGrm {
    executor 'slurm'
    queue 'gr10478b'
    time '48h'
    publishDir "${params.resultsDir}/02.sparseGrm/02.nullModel", mode: 'symlink'

    input:
    tuple path(sparse_grm), path(sparse_grm_id)
    tuple val(plink_prefix), path(plink_bed), path(plink_bim), path(plink_fam)
    path(pheno_file)
    path(cov_file)
    path(age_na_remove_file)

    output:
    tuple path('*.rda'), path('*.varianceRatio.txt'), emit: nullModel

    script:
    def mergeScript = "${params.scriptDir}/merge_pheno_cov.py"
    def covarTagScript = "${params.scriptDir}/covar_tag_from_names.py"
    def plinkPrefix = plink_bed.baseName
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.condaEnvActivate}
    set -euo pipefail
    covarTag=\$(python3 ${covarTagScript} --covar-names "${params.covarNames}")
    sparseNullPrefix="${params.projectName}.\${covarTag}.sparseGrm.null"

    # Merge phenotype and covariates once for Step 1
    python3 ${mergeScript} \
        --pheno "${pheno_file}" \\
        --cov "${cov_file}" \\
        --pheno_col "${params.phenoName}" \\
        --cov_list "${params.covarNames}" \
        --sex_col "SEX" \\
        --out "merged_pheno_cov.txt"

    if [[ "${params.includeAgeZ}" == "true" ]]; then
        awk 'NR==FNR{rm[\$2]=1; next} FNR==1 || !(\$2 in rm)' ${age_na_remove_file} merged_pheno_cov.txt > merged_pheno_cov.filtered.txt
    else
        cp merged_pheno_cov.txt merged_pheno_cov.filtered.txt
    fi

    # SAIGE Step 1 (sparse GRM)
    step1_fitNULLGLMM.R \
        --plinkFile=${plinkPrefix} \
        --phenoFile=merged_pheno_cov.filtered.txt \
        --phenoCol=${params.phenoName} \
        --covarColList=${params.covarNames} \
        --sexCol=SEX \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --outputPrefix=\${sparseNullPrefix} \
        --nThreads=32 \
        --numRandomMarkerforVarianceRatio=200 \
        --skipVarianceRatioEstimation=FALSE \
        --useSparseGRMtoFitNULL=TRUE  \
        --sparseGRMFile=${sparse_grm} \
        --sparseGRMSampleIDFile=${sparse_grm_id} \
        --IsOverwriteVarianceRatioFile=TRUE \
        --isCovariateOffset=FALSE
    """
}

/*
 * Stage 3 (Sparse): Run per-chromosome association tests.
 */
process AssocTestSparseGrm {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    publishDir "${params.resultsDir}/02.sparseGrm/03.assocByChr", mode: 'symlink'

    input:
    tuple path(model_file), path(variance_ratio), path(sparse_grm_file), path(sparse_grm_id)
    tuple val(plink_prefix), path(plink_bed), path(plink_bim), path(plink_fam)
    each chr

    output:
    tuple val(chr), path('*.assoc.txt'), emit: assocResults

    script:
    def covarTagScript = "${params.scriptDir}/covar_tag_from_names.py"

    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.condaEnvActivate}
    set -euo pipefail
    covarTag=\$(python3 ${covarTagScript} --covar-names "${params.covarNames}")
    outFile="${params.projectName}.\${covarTag}.sparseGrm.chr${chr}.assoc.txt"

    # SAIGE Step 2 (sparse GRM)
    step2_SPAtests.R \
        --bedFile=${plink_bed} \
        --bimFile=${plink_bim} \
        --famFile=${plink_fam} \
        --chrom=${chr} \
        --AlleleOrder=alt-first \
        --SAIGEOutputFile=\$outFile \
        --GMMATmodelFile=${model_file} \
        --varianceRatioFile=${variance_ratio} \
        --LOCO=FALSE \
        --SPAcutoff=1.645 \
        --is_output_moreDetails=TRUE \
        --sparseGRMFile=${sparse_grm_file} \
        --sparseGRMSampleIDFile=${sparse_grm_id} \
        --is_fastTest=FALSE
    """
}

/*
 * Stage 4 (Shared): Merge chr-wise association outputs.
 */
process MergeAssocResults {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    publishDir "${params.resultsDir}/${analysisName}/04.assocMerged", mode: 'symlink'

    input:
    val(analysisName)
    path(assoc_files)

    output:
    tuple val(analysisName), path('*.assoc.merged.txt'), emit: mergedResults

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.condaEnvActivate}
    set -euo pipefail

    covarTag=\$(python3 ${params.scriptDir}/covar_tag_from_names.py --covar-names "${params.covarNames}")
    if [[ "${analysisName}" == *"full"* ]]; then
        grmType="fullGrm"
    else
        grmType="sparseGrm"
    fi
    mergedOut="${params.projectName}.\${covarTag}.\${grmType}.assoc.merged.txt"

    assocCount=\$(ls -1 *.assoc.txt 2>/dev/null | wc -l)
    if [[ "\${assocCount}" -eq 0 ]]; then
        echo "Error: No per-chromosome association files found for \${analysisName}" >&2
        exit 1
    fi

    firstFile=\$(ls *.assoc.txt | sort -V | head -n 1)
    head -n 1 "\${firstFile}" > \${mergedOut}
    for file in \$(ls *.assoc.txt | sort -V); do
        tail -n +2 "\$file" >> \${mergedOut}
    done
    """
}

/*
 * Stage 5 (Shared): Generate Manhattan and QQ plots from merged results.
 */
process PlotManhattanQq {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    publishDir "${params.resultsDir}/${analysisName}/05.plots", mode: 'symlink'

    input:
    tuple val(analysisName), path(assocFile)
    val(plotTitle)

    output:
    tuple val(analysisName), path('*.png'), emit: plots

    script:
    def scriptPath = "${params.scriptDir}/saige_manhattan_qq.py"
    def outPrefix = assocFile.baseName
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.condaEnvActivate}
    set -euo pipefail

    python3 ${scriptPath} \
        --input ${assocFile} \
        --output-prefix ${outPrefix} \
        --title "${plotTitle}"
    """
}

/*
 * Post-processing subworkflow for Full-GRM branch.
 * Using a separate workflow context avoids DSL2 component reuse conflicts.
 */
workflow FullGrmPost {
    take:
    assocFiles

    main:
    MergeAssocResults('01.fullGrm', assocFiles)
    PlotManhattanQq(
        MergeAssocResults.out.mergedResults,
        channel.value('SAIGE Results (Full GRM)')
    )
}

/*
 * Post-processing subworkflow for Sparse-GRM branch.
 * Using a separate workflow context avoids DSL2 component reuse conflicts.
 */
workflow SparseGrmPost {
    take:
    assocFiles

    main:
    MergeAssocResults('02.sparseGrm', assocFiles)
    PlotManhattanQq(
        MergeAssocResults.out.mergedResults,
        channel.value('SAIGE Results (Sparse GRM)')
    )
}

/*
 * Main Workflow
 */
workflow {
    // Assemble canonical input paths.
    plinkPrefix = "${params.genotypeRoot}/${params.randomModelDirRelPath}/${params.plinkPrefixName}"
    phenoPath = "${params.genotypeRoot}/${params.covPhenoDirRelPath}/${params.phenoFileName}"
    covPath = "${params.genotypeRoot}/${params.covPhenoDirRelPath}/${params.covFileName}"
    ageNaRemovePath = "${params.genotypeRoot}/${params.covPhenoDirRelPath}/${params.ageNaRemoveFileName}"
    highLdPath = "${params.infoDir}/${params.highLdFileName}"

    // Early validation for key parameters and required input files.
    if (!(params.covarPcSource in ['bbj', 'own'])) {
        error "Invalid params.covarPcSource='${params.covarPcSource}'. Allowed values: bbj | own"
    }
    [
        "${plinkPrefix}.bed",
        "${plinkPrefix}.bim",
        "${plinkPrefix}.fam",
        phenoPath,
        covPath,
        highLdPath,
        "${params.scriptDir}/normalize_high_ld_chr_prefix.sh",
        "${params.scriptDir}/covar_tag_from_names.py",
        "${params.scriptDir}/merge_pheno_cov.py",
        "${params.scriptDir}/saige_manhattan_qq.py",
    ].each { reqPath ->
        if (!file(reqPath).exists()) {
            error "Missing required input/script path: ${reqPath}"
        }
    }
    if (params.includeAgeZ && !file(ageNaRemovePath).exists()) {
        error "Missing AGE_Z removal file: ${ageNaRemovePath}"
    }
    
    // Input channels
    plinkCh = channel.fromFilePairs("${plinkPrefix}.{bed,bim,fam}", size: 3)
        .map { id, files -> tuple(id, files[0], files[1], files[2]) } 
    
    // Channels for Phenotype and Covariate files
    phenoCh = channel.fromPath(phenoPath)
    covCh = channel.fromPath(covPath)
    ageNaRemoveCh = channel.fromPath(ageNaRemovePath)

    // Autosomes
    chrCh = channel.of(1..22)
    
    // High LD regions file
    highLdCh = channel.fromPath(highLdPath)

    // Shared stage
    targetPlinkCh = plinkCh
    if (params.includeAgeZ) {
        RemoveAgeNaSamples(plinkCh, ageNaRemoveCh)
        targetPlinkCh = RemoveAgeNaSamples.out.filteredPlink
    }

    LdPruning(targetPlinkCh, highLdCh)
    pruneInCh = LdPruning.out.pruneIn

    // Full-GRM branch (mandatory)
    SplitPlinkToBgen(targetPlinkCh, chrCh)
    FitNullGlmmFullGrm(targetPlinkCh, pruneInCh, phenoCh, covCh, ageNaRemoveCh)

    fullAssocInputCh = FitNullGlmmFullGrm.out.nullModel.combine(SplitPlinkToBgen.out.bgenFiles)
    AssocTestFullGrm(fullAssocInputCh)
    fullAssocFilesCh = AssocTestFullGrm.out.assocResults
        .map { row -> row[1] }
        .collect()
    FullGrmPost(fullAssocFilesCh)

    // Sparse-GRM branch (optional)
    if (params.runSparse) {
        CalcSparseGrm(targetPlinkCh, pruneInCh)
        FitNullGlmmSparseGrm(
            CalcSparseGrm.out.sparseGrm,
            targetPlinkCh,
            phenoCh,
            covCh,
            ageNaRemoveCh
        )

        sparseAssocInputCh = FitNullGlmmSparseGrm.out.nullModel.combine(CalcSparseGrm.out.sparseGrm)
        AssocTestSparseGrm(
            sparseAssocInputCh,
            targetPlinkCh,
            chrCh
        )

        sparseAssocFilesCh = AssocTestSparseGrm.out.assocResults
            .map { row -> row[1] }
            .collect()
        SparseGrmPost(sparseAssocFilesCh)
    }
}