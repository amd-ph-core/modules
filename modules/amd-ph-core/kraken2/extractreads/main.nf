// KRAKEN2_EXTRACTREADS: Kraken2 classification and read extraction in one task.
//
// Running both steps in one task means the per-read assignment file (one line per
// read) is read where it is written, instead of being staged into a second task.
// Extraction uses extract_kraken_reads_parallel.py, which ships in this module's
// resources/usr/bin; pipelines need `nextflow.enable.moduleBinaries = true`.

process KRAKEN2_EXTRACTREADS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f827dcea51be6b5c32255167caa2dfb65607caecdc8b067abd6b71c267e2e82/data' :
        'community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71' }"

    input:
    tuple val(meta), path(reads)
    path db
    val taxid
    val save_output_fastqs
    val save_reads_assignment

    output:
    tuple val(meta), path('*.extracted_kraken2_read*.fastq.gz')  , emit: extracted_kraken2_reads
    tuple val(meta), path('*.kraken2.report.txt')                 , emit: report
    tuple val(meta), env('READS_PRESENT')                         , emit: reads_present
    tuple val(meta), path('*.extract_kraken_reads.status.json')   , emit: status
    tuple val(meta), path('*.classified{.,_}*.fastq.gz')          , emit: classified_reads_fastq     , optional: true
    tuple val(meta), path('*.unclassified{.,_}*.fastq.gz')        , emit: unclassified_reads_fastq   , optional: true
    tuple val(meta), path('*.kraken2.classifiedreads.txt')        , emit: classified_reads_assignment, optional: true
    tuple val("${task.process}"), val('kraken2'), eval('kraken2 --version 2>&1 | head -1 | sed "s/^.*Kraken version //; s/ .*//"'), emit: versions_kraken2, topic: versions
    tuple val("${task.process}"), val('extract_kraken_reads_parallel.py'), eval('extract_kraken_reads_parallel.py --version'), emit: versions_extractor, topic: versions
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/^pigz //"'), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args  ?: ''   // kraken2
    def args2  = task.ext.args2 ?: ''   // extract_kraken_reads_parallel.py, e.g. --exclude
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? '' : '--paired'
    def classified   = meta.single_end ? "${prefix}.classified.fastq"   : "${prefix}.classified#.fastq"
    def unclassified = meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    def classified_option   = save_output_fastqs ? "--classified-out ${classified}" : ''
    def unclassified_option = save_output_fastqs ? "--unclassified-out ${unclassified}" : ''
    def read_list    = reads instanceof List ? reads : [reads]
    def input_reads  = meta.single_end ? "-s ${read_list[0]}" : "-s1 ${read_list[0]} -s2 ${read_list[1]}"
    def output_reads = meta.single_end ?
        "-o ${prefix}.extracted_kraken2_read.fastq.gz" :
        "-o ${prefix}.extracted_kraken2_read_1.fastq.gz -o2 ${prefix}.extracted_kraken2_read_2.fastq.gz"
    def kraken_threads = Math.max(1, task.cpus - 1)
    """
    kraken2 \\
        --db ${db} \\
        --threads ${kraken_threads} \\
        --report ${prefix}.kraken2.report.txt \\
        --gzip-compressed \\
        --output ${prefix}.kraken2.classifiedreads.txt \\
        ${unclassified_option} \\
        ${classified_option} \\
        ${paired} \\
        ${args} \\
        ${reads}

    extract_kraken_reads_parallel.py \\
        ${args2} \\
        --threads ${task.cpus} \\
        --gzip-output \\
        --status-file ${prefix}.extract_kraken_reads.status.json \\
        --fastq-output \\
        -t ${taxid} \\
        -k ${prefix}.kraken2.classifiedreads.txt \\
        -r ${prefix}.kraken2.report.txt \\
        ${input_reads} \\
        ${output_reads}

    READS_PRESENT=\$(python3 -c "import json; print(str(json.load(open('${prefix}.extract_kraken_reads.status.json'))['reads_present']).lower())")

    if [ "${save_output_fastqs}" = "true" ]; then
        pigz -p ${task.cpus} -- *.fastq
    fi
    if [ "${save_reads_assignment}" != "true" ]; then
        rm -f ${prefix}.kraken2.classifiedreads.txt
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extracted = meta.single_end ?
        "${prefix}.extracted_kraken2_read.fastq.gz" :
        "${prefix}.extracted_kraken2_read_1.fastq.gz ${prefix}.extracted_kraken2_read_2.fastq.gz"
    def classified = meta.single_end ?
        "${prefix}.classified.fastq.gz ${prefix}.unclassified.fastq.gz" :
        "${prefix}.classified_1.fastq.gz ${prefix}.classified_2.fastq.gz ${prefix}.unclassified_1.fastq.gz ${prefix}.unclassified_2.fastq.gz"
    """
    touch ${prefix}.kraken2.report.txt
    for output in ${extracted}; do
        printf '' | gzip -n > "\$output"
    done
    if [ "${save_output_fastqs}" = "true" ]; then
        for output in ${classified}; do
            printf '' | gzip -n > "\$output"
        done
    fi
    if [ "${save_reads_assignment}" = "true" ]; then
        touch ${prefix}.kraken2.classifiedreads.txt
    fi
    echo '{"reads_present": true}' > ${prefix}.extract_kraken_reads.status.json
    READS_PRESENT=true
    """
}
