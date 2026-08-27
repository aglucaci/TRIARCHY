process WGS_QC {
    tag "${library_id}"
    label 'medium'

    input:
    tuple val(library_id), val(layout), val(r1), val(r2)

    output:
    tuple val(library_id), val(layout), path("${library_id}"), emit: reads
    tuple val(library_id), path("${library_id}.fastp.json"), emit: metrics

    script:
    def pairedArgs = layout == 'paired' ?
        "-i '${r1}' -I '${r2}' -o '${library_id}/R1.fastq.gz' -O '${library_id}/R2.fastq.gz'" :
        "-i '${r1}' -o '${library_id}/R1.fastq.gz'"
    """
    command -v fastp >/dev/null || { echo 'fastp is required' >&2; exit 127; }
    mkdir -p '${library_id}'
    fastp ${pairedArgs} \
        --json '${library_id}/fastp.json' \
        --html '${library_id}/fastp.html' \
        --thread ${task.cpus} \
        ${params.fastp_extra}
    cp '${library_id}/fastp.json' '${library_id}.fastp.json'
    """

}

process WGS_HOST_REMOVAL {
    tag "${library_id}"
    label 'high'

    input:
    tuple val(library_id), val(layout), path(qc_dir)

    output:
    tuple val(library_id), val(layout), path("${library_id}"), emit: reads
    tuple val(library_id), path("${library_id}.bowtie2.log"), emit: metrics

    script:
    if (!params.host_index) {
        error "Production WMS execution requires --host_index"
    }
    def hostIndex = file(params.host_index, checkIfExists: false)
    def readArgs = layout == 'paired' ?
        "-1 '${qc_dir}/R1.fastq.gz' -2 '${qc_dir}/R2.fastq.gz' --un-conc-gz '${library_id}/clean.fastq.gz'" :
        "-U '${qc_dir}/R1.fastq.gz' --un-gz '${library_id}/R1.fastq.gz'"
    def rename = layout == 'paired' ?
        "mv '${library_id}/clean.fastq.1.gz' '${library_id}/R1.fastq.gz'; mv '${library_id}/clean.fastq.2.gz' '${library_id}/R2.fastq.gz'" :
        "true"
    """
    command -v bowtie2 >/dev/null || { echo 'bowtie2 is required' >&2; exit 127; }
    mkdir -p '${library_id}'
    bowtie2 -x '${hostIndex}' ${readArgs} \
        --threads ${task.cpus} \
        --met-file '${library_id}/bowtie2.metrics.txt' \
        -S /dev/null 2> '${library_id}/bowtie2.log'
    cp '${library_id}/bowtie2.log' '${library_id}.bowtie2.log'
    ${rename}
    """

}

process WGS_TAXONOMY {
    tag "${library_id}"
    label 'high'

    input:
    tuple val(library_id), val(layout), path(clean_dir)

    output:
    tuple val(library_id), path("${library_id}.taxonomy/profile.tsv"), emit: profile
    path "${library_id}.taxonomy", emit: native_files

    script:
    if (!params.kraken2_db) {
        error "Production WMS execution requires --kraken2_db"
    }
    def krakenDb = file(params.kraken2_db, checkIfExists: true)
    def reads = layout == 'paired' ?
        "--paired '${clean_dir}/R1.fastq.gz' '${clean_dir}/R2.fastq.gz'" :
        "'${clean_dir}/R1.fastq.gz'"
    def memoryMappingArg = params.kraken2_memory_mapping ? '--memory-mapping' : ''
    """
    command -v kraken2 >/dev/null || { echo 'kraken2 is required' >&2; exit 127; }
    command -v bracken >/dev/null || { echo 'bracken is required' >&2; exit 127; }
    mkdir -p '${library_id}.taxonomy'
    kraken2 --db '${krakenDb}' ${memoryMappingArg} ${reads} \
        --threads ${task.cpus} \
        --report '${library_id}.taxonomy/kraken2.report' \
        --output '${library_id}.taxonomy/kraken2.output'
    bracken -d '${krakenDb}' \
        -i '${library_id}.taxonomy/kraken2.report' \
        -o '${library_id}.taxonomy/bracken.species.tsv' \
        -r ${params.bracken_read_length} -l S
    python3 ${projectDir}/workflow/scripts/triarchy_integrate.py standardize-bracken \
        --library-id '${library_id}' \
        --input '${library_id}.taxonomy/bracken.species.tsv' \
        --output '${library_id}.taxonomy/profile.tsv'
    """

}

process WGS_FUNCTIONS {
    tag "${library_id}"
    label 'high'

    input:
    tuple val(library_id), val(layout), path(clean_dir)

    output:
    tuple val(library_id), path("${library_id}.functions/functions.tsv"), emit: profile
    path "${library_id}.functions", emit: native_files

    script:
    def prep = layout == 'paired' ?
        "cat '${clean_dir}/R1.fastq.gz' '${clean_dir}/R2.fastq.gz' > '${library_id}.functions/combined.fastq.gz'" :
        "cp '${clean_dir}/R1.fastq.gz' '${library_id}.functions/combined.fastq.gz'"
    def humannDb = params.humann_db ?
        file(params.humann_db, checkIfExists: true) : null
    def humannProteinDb = params.humann_protein_db ?
        file(params.humann_protein_db, checkIfExists: true) : null
    def metaphlanDb = params.metaphlan_db ?
        file(params.metaphlan_db, checkIfExists: true) : null
    if (!params.humann_bypass_prescreen && (!metaphlanDb || !params.metaphlan_index)) {
        error "Production HUMAnN execution requires --metaphlan_db and --metaphlan_index"
    }
    def nucleotideArg = humannDb ? "--nucleotide-database '${humannDb}'" : ''
    def proteinArg = humannProteinDb ? "--protein-database '${humannProteinDb}'" : ''
    def prescreenArg = params.humann_bypass_prescreen ? '--bypass-prescreen' : ''
    def metaphlanArg = metaphlanDb ?
        "--metaphlan-options '-t rel_ab --db_dir ${metaphlanDb} -x ${params.metaphlan_index} --offline'" : ''
    """
    command -v humann >/dev/null || { echo 'humann is required' >&2; exit 127; }
    metaphlan_real="\$(command -v metaphlan)" || { echo 'metaphlan is required by HUMAnN' >&2; exit 127; }
    command -v diamond >/dev/null || { echo 'diamond is required by HUMAnN' >&2; exit 127; }
    diamond version >/dev/null || { echo 'diamond is installed but cannot execute on this platform' >&2; exit 126; }
    export TRIARCHY_METAPHLAN_REAL="\${metaphlan_real}"
    export PATH="${projectDir}/workflow/bin:\${PATH}"
    mkdir -p '${library_id}.functions'
    ${prep}
    humann --input '${library_id}.functions/combined.fastq.gz' \
        --output '${library_id}.functions/humann' \
        --threads ${task.cpus} ${nucleotideArg} ${proteinArg} ${prescreenArg} ${metaphlanArg}
    python3 ${projectDir}/workflow/scripts/triarchy_integrate.py standardize-humann \
        --library-id '${library_id}' \
        --input-dir '${library_id}.functions/humann' \
        --output '${library_id}.functions/functions.tsv'
    """

}
