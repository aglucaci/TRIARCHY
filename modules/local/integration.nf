process QC_SUMMARY {
    tag 'assay quality'
    label 'low'
    publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

    input:
    path samples
    path qc_files

    output:
    path 'library_qc.tsv', emit: metrics

    script:
    """
    python3 ${projectDir}/workflow/scripts/triarchy_integrate.py summarize-qc \
        --samples '${samples}' \
        --qc-dirs ${qc_files.join(' ')} \
        --output library_qc.tsv
    """
}

process CONTROL_ANNOTATION {
    tag 'linked controls'
    label 'low'
    publishDir "${params.outdir}/profiles", mode: 'copy', overwrite: true

    input:
    path samples
    path controls
    path taxonomic_dirs, stageAs: 'taxa/*'
    path function_dirs, stageAs: 'functions/*'

    output:
    path 'unfiltered/taxa.tsv', emit: unfiltered_taxa
    path 'unfiltered/functions.tsv', emit: unfiltered_functions
    path 'control_annotated/taxa.tsv', emit: annotated_taxa
    path 'control_annotated/functions.tsv', emit: annotated_functions
    path 'high_confidence/taxa.tsv', emit: high_confidence
    path 'sample_confidence.tsv', emit: sample_confidence

    script:
    """
    mkdir -p unfiltered control_annotated high_confidence
    python3 ${projectDir}/workflow/scripts/triarchy_integrate.py control-annotate \
        --samples '${samples}' \
        --controls '${controls}' \
        --taxa ${taxonomic_dirs.join(' ')} \
        --functions ${function_dirs.join(' ')} \
        --unfiltered-taxa unfiltered/taxa.tsv \
        --unfiltered-functions unfiltered/functions.tsv \
        --annotated-taxa control_annotated/taxa.tsv \
        --annotated-functions control_annotated/functions.tsv \
        --high-confidence high_confidence/taxa.tsv \
        --sample-confidence sample_confidence.tsv
    """
}

process TAXONOMY_HARMONIZATION {
    tag 'taxonomy'
    label 'low'
    publishDir "${params.outdir}/integrated", mode: 'copy', overwrite: true

    input:
    path annotated_taxa
    val crosswalk

    output:
    path 'taxonomy_crosswalk.tsv', emit: crosswalk
    path 'harmonized_taxa.tsv', emit: harmonized

    script:
    def crosswalkArg = crosswalk ? "--crosswalk '${crosswalk}'" : ''
    """
    python3 ${projectDir}/workflow/scripts/triarchy_integrate.py harmonize \
        --input '${annotated_taxa}' \
        ${crosswalkArg} \
        --output harmonized_taxa.tsv \
        --output-crosswalk taxonomy_crosswalk.tsv
    """
}

process EVIDENCE_CLASSIFICATION {
    tag 'cross-assay evidence'
    label 'low'
    publishDir "${params.outdir}/integrated", mode: 'copy', overwrite: true

    input:
    path harmonized
    path samples
    path assay_performance
    path evidence_policy
    path qc_metrics

    output:
    path 'taxon_evidence.tsv', emit: evidence
    path 'evidence_model.json', emit: model

    script:
    """
    python3 ${projectDir}/workflow/scripts/triarchy_integrate.py classify \
        --harmonized '${harmonized}' \
        --samples '${samples}' \
        --assay-performance '${assay_performance}' \
        --evidence-policy '${evidence_policy}' \
        --qc-metrics '${qc_metrics}' \
        --model-output evidence_model.json \
        --output taxon_evidence.tsv
    """
}

process CROSS_ASSAY_CONCORDANCE {
    tag 'cross-assay concordance'
    label 'low'
    publishDir "${params.outdir}/integrated", mode: 'copy', overwrite: true

    input:
    path harmonized
    path samples

    output:
    path 'cross_assay_concordance.tsv', emit: concordance

    script:
    """
    python3 ${projectDir}/workflow/scripts/triarchy_integrate.py concordance \
        --harmonized '${harmonized}' \
        --samples '${samples}' \
        --output cross_assay_concordance.tsv
    """
}

process ASSAY_DIVERSITY {
    tag 'assay-specific diversity'
    label 'low'
    publishDir "${params.outdir}/diversity", mode: 'copy', overwrite: true

    input:
    path high_confidence
    path annotated_functions

    output:
    path 'alpha_diversity.tsv', emit: alpha
    path 'beta_diversity.tsv', emit: beta

    script:
    """
    python3 ${projectDir}/workflow/scripts/triarchy_integrate.py diversity \
        --taxa '${high_confidence}' \
        --functions '${annotated_functions}' \
        --alpha-output alpha_diversity.tsv \
        --beta-output beta_diversity.tsv
    """
}

process RUN_PROVENANCE {
    tag 'provenance'
    label 'low'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: true

    input:
    path samples
    path controls
    path validated

    output:
    path 'provenance.json', emit: provenance

    script:
    """
    python3 ${projectDir}/workflow/scripts/triarchy_integrate.py provenance \
        --samples '${samples}' \
        --controls '${controls}' \
        --project-dir '${projectDir}' \
        --output provenance.json
    """
}

process BUILD_REPORT {
    tag 'HTML report'
    label 'low'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path evidence
    path concordance
    path alpha_diversity
    path beta_diversity
    path sample_confidence
    path evidence_model
    path provenance

    output:
    path 'triarchy_report.html', emit: report

    script:
    """
    python3 ${projectDir}/workflow/scripts/triarchy_integrate.py report \
        --evidence '${evidence}' \
        --concordance '${concordance}' \
        --alpha-diversity '${alpha_diversity}' \
        --beta-diversity '${beta_diversity}' \
        --sample-confidence '${sample_confidence}' \
        --evidence-model '${evidence_model}' \
        --provenance '${provenance}' \
        --output triarchy_report.html
    """
}
