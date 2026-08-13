process VALIDATE_METADATA {
    tag 'metadata'
    label 'low'

    input:
    path samples
    path controls

    output:
    path 'metadata.validated', emit: validated

    script:
    """
    python3 ${projectDir}/scripts/validate_benchmark_metadata.py \
        '${samples}' '${controls}'
    touch metadata.validated
    """
}

process SPLIT_MANIFESTS {
    tag 'assay manifests'
    label 'low'

    input:
    path samples
    path validated
    val samples_base_dir

    output:
    path 'wms.tsv', emit: wms
    path '16s.tsv', emit: s16
    path 'its.tsv', emit: its

    script:
    """
    python3 ${projectDir}/workflow/scripts/triarchy_integrate.py split-manifest \
        --samples '${samples}' \
        --base-dir '${samples_base_dir}' \
        --outdir .
    """
}
