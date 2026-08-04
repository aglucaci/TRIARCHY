nextflow.enable.dsl = 2

include { VALIDATE_METADATA; SPLIT_MANIFESTS } from './modules/local/metadata'
include { WMS_QC; WMS_HOST_REMOVAL; WMS_TAXONOMY; WMS_FUNCTIONS } from './modules/local/wms'
include { AMPLICON_QC as S16_QC; DADA2_PROFILE as S16_DADA2 } from './modules/local/amplicon'
include { AMPLICON_QC as ITS_QC; ITSXPRESS; DADA2_PROFILE as ITS_DADA2 } from './modules/local/amplicon'
include {
    QC_SUMMARY
    CONTROL_ANNOTATION
    TAXONOMY_HARMONIZATION
    EVIDENCE_CLASSIFICATION
    CROSS_ASSAY_CONCORDANCE
    ASSAY_DIVERSITY
    RUN_PROVENANCE
    BUILD_REPORT
} from './modules/local/integration'

workflow {
    if (!params.samples || !params.controls) {
        error 'Required parameters: --samples and --controls'
    }
    samples = file(params.samples, checkIfExists: true)
    controls = file(params.controls, checkIfExists: true)
    s16_taxonomy = params.s16_taxonomy ?
        file(params.s16_taxonomy, checkIfExists: true) : ''
    its_taxonomy = params.its_taxonomy ?
        file(params.its_taxonomy, checkIfExists: true) : ''
    taxonomy_crosswalk = params.taxonomy_crosswalk ?
        file(params.taxonomy_crosswalk, checkIfExists: true) : ''
    assay_performance = file(
        params.assay_performance ?: "${projectDir}/conf/assay_performance.tsv",
        checkIfExists: true
    )
    evidence_policy = file(
        params.evidence_policy ?: "${projectDir}/conf/evidence_policy.json",
        checkIfExists: true
    )

    VALIDATE_METADATA(samples, controls)
    SPLIT_MANIFESTS(
        samples,
        VALIDATE_METADATA.out.validated,
        samples.parent.toString()
    )

    wms_rows = SPLIT_MANIFESTS.out.wms
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.library_id,
                row.library_layout,
                row.fastq_r1,
                row.fastq_r2 ?: ''
            )
        }
    s16_rows = SPLIT_MANIFESTS.out.s16
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.library_id,
                row.library_layout,
                row.fastq_r1,
                row.fastq_r2 ?: ''
            )
        }
    its_rows = SPLIT_MANIFESTS.out.its
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.library_id,
                row.library_layout,
                row.fastq_r1,
                row.fastq_r2 ?: ''
            )
        }

    WMS_QC(wms_rows)
    WMS_HOST_REMOVAL(WMS_QC.out.reads)
    WMS_TAXONOMY(WMS_HOST_REMOVAL.out.reads)
    WMS_FUNCTIONS(WMS_HOST_REMOVAL.out.reads)

    S16_QC(s16_rows, '16S')
    S16_DADA2(S16_QC.out.reads, s16_taxonomy)

    ITS_QC(its_rows, 'ITS')
    ITSXPRESS(ITS_QC.out.reads)
    ITS_DADA2(ITSXPRESS.out.reads, its_taxonomy)

    qc_files = WMS_QC.out.metrics
        .mix(WMS_HOST_REMOVAL.out.metrics)
        .mix(S16_QC.out.metrics)
        .mix(ITS_QC.out.metrics)
        .map { library_id, metric -> metric }
        .collect()
    QC_SUMMARY(samples, qc_files)

    taxonomic_profiles = WMS_TAXONOMY.out.native_files
        .mix(S16_DADA2.out.native_files)
        .mix(ITS_DADA2.out.native_files)
        .collect()
    function_profiles = WMS_FUNCTIONS.out.native_files.collect()

    CONTROL_ANNOTATION(
        samples,
        controls,
        taxonomic_profiles,
        function_profiles
    )
    TAXONOMY_HARMONIZATION(
        CONTROL_ANNOTATION.out.annotated_taxa,
        taxonomy_crosswalk
    )
    EVIDENCE_CLASSIFICATION(
        TAXONOMY_HARMONIZATION.out.harmonized,
        samples,
        assay_performance,
        evidence_policy,
        QC_SUMMARY.out.metrics
    )
    CROSS_ASSAY_CONCORDANCE(
        TAXONOMY_HARMONIZATION.out.harmonized,
        samples
    )
    ASSAY_DIVERSITY(
        CONTROL_ANNOTATION.out.high_confidence,
        CONTROL_ANNOTATION.out.annotated_functions
    )
    RUN_PROVENANCE(
        samples,
        controls,
        VALIDATE_METADATA.out.validated
    )
    BUILD_REPORT(
        EVIDENCE_CLASSIFICATION.out.evidence,
        CROSS_ASSAY_CONCORDANCE.out.concordance,
        ASSAY_DIVERSITY.out.alpha,
        ASSAY_DIVERSITY.out.beta,
        CONTROL_ANNOTATION.out.sample_confidence,
        EVIDENCE_CLASSIFICATION.out.model,
        RUN_PROVENANCE.out.provenance
    )
}
