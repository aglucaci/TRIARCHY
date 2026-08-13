process AMPLICON_QC {
    tag "${assay}:${library_id}"
    label 'medium'

    input:
    tuple val(library_id), val(layout), val(r1), val(r2)
    val assay

    output:
    tuple val(library_id), val(layout), val(assay), path("${library_id}"), emit: reads
    tuple val(library_id), path("${library_id}.cutadapt.json"), emit: metrics

    script:
    def forward = assay == '16S' ? params.cutadapt_16s_forward : params.cutadapt_its_forward
    def reverse = assay == '16S' ? params.cutadapt_16s_reverse : params.cutadapt_its_reverse
    if (!forward) {
        error "Production ${assay} execution requires primer parameters in nextflow.config or on the command line"
    }
    def paired = layout == 'paired'
    def inputArgs = paired ? "'${r1}' '${r2}'" : "'${r1}'"
    def outputArgs = paired ?
        "-o '${library_id}/R1.fastq.gz' -p '${library_id}/R2.fastq.gz'" :
        "-o '${library_id}/R1.fastq.gz'"
    def adapterArgs = paired && reverse ?
        "-g '${forward}' -G '${reverse}'" :
        "-g '${forward}'"
    """
    command -v cutadapt >/dev/null || { echo 'cutadapt is required' >&2; exit 127; }
    mkdir -p '${library_id}'
    cutadapt ${adapterArgs} ${outputArgs} ${inputArgs} \
        --cores ${task.cpus} \
        --json '${library_id}/cutadapt.json'
    cp '${library_id}/cutadapt.json' '${library_id}.cutadapt.json'
    """

    stub:
    """
    mkdir -p '${library_id}'
    touch '${library_id}/R1.fastq.gz'
    if [[ '${layout}' == 'paired' ]]; then touch '${library_id}/R2.fastq.gz'; fi
    printf '{"stub":true}\n' > '${library_id}/cutadapt.json'
    cp '${library_id}/cutadapt.json' '${library_id}.cutadapt.json'
    """
}

process ITSXPRESS {
    tag "${library_id}"
    label 'medium'

    input:
    tuple val(library_id), val(layout), val(assay), path(qc_dir)

    output:
    tuple val(library_id), val(layout), val(assay), path("${library_id}.itsxpress"), emit: reads

    script:
    def pairedArgs = layout == 'paired' ?
        "--fastq '${qc_dir}/R1.fastq.gz' --fastq2 '${qc_dir}/R2.fastq.gz' --outfile '${library_id}.itsxpress/R1.fastq.gz' --outfile2 '${library_id}.itsxpress/R2.fastq.gz'" :
        "--fastq '${qc_dir}/R1.fastq.gz' --outfile '${library_id}.itsxpress/R1.fastq.gz'"
    def emptyCopy = layout == 'paired' ?
        "cp '${qc_dir}/R1.fastq.gz' '${library_id}.itsxpress/R1.fastq.gz'; cp '${qc_dir}/R2.fastq.gz' '${library_id}.itsxpress/R2.fastq.gz'" :
        "cp '${qc_dir}/R1.fastq.gz' '${library_id}.itsxpress/R1.fastq.gz'"
    """
    command -v itsxpress >/dev/null || { echo 'itsxpress is required' >&2; exit 127; }
    mkdir -p '${library_id}.itsxpress'
    if ! python3 -c "import gzip, sys; handle=gzip.open(sys.argv[1], 'rb'); sys.exit(0 if handle.read(1) else 1)" '${qc_dir}/R1.fastq.gz'; then
        ${emptyCopy}
        exit 0
    fi
    itsxpress ${pairedArgs} --region ITS2 --taxa All --threads ${task.cpus}
    """

    stub:
    """
    mkdir -p '${library_id}.itsxpress'
    touch '${library_id}.itsxpress/R1.fastq.gz'
    if [[ '${layout}' == 'paired' ]]; then touch '${library_id}.itsxpress/R2.fastq.gz'; fi
    """
}

process DADA2_PROFILE {
    tag "${assay}:${library_id}"
    label 'high'

    input:
    tuple val(library_id), val(layout), val(assay), path(read_dir)
    val taxonomy_reference

    output:
    tuple val(library_id), path("${library_id}/profile.tsv"), emit: profile
    path "${library_id}", emit: native_files

    script:
    if (!taxonomy_reference) {
        error "Production ${assay} execution requires a taxonomy reference"
    }
    """
    command -v Rscript >/dev/null || { echo 'Rscript is required' >&2; exit 127; }
    mkdir -p '${library_id}'
    Rscript ${projectDir}/workflow/scripts/dada2_profile.R \
        '${library_id}' '${assay}' '${layout}' \
        '${read_dir}/R1.fastq.gz' '${read_dir}/R2.fastq.gz' \
        '${taxonomy_reference}' '${library_id}'
    """

    stub:
    """
    mkdir -p '${library_id}'
    python3 ${projectDir}/workflow/scripts/triarchy_integrate.py mock-profile \
        --library-id '${library_id}' --assay '${assay}' \
        --output '${library_id}/profile.tsv'
    printf '>ASV1\\nACGTACGT\\n' > '${library_id}/asvs.fasta'
    printf 'feature_id\\tcount\\nASV1\\t100\\n' > '${library_id}/asv_counts.tsv'
    """
}
