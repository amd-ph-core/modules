process PICARD_SORTSAM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'cdc-amd/picard:3.4.0_ac49616_v8'

    input:
    tuple val(meta), path(bam)
    val sort_order

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('picard'), eval('picard SortSam --version 2>&1 | grep -o "Version:.*" | cut -f2- -d:'), emit: versions_picard, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard SortSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    picard \\
        SortSam \\
        -Xmx${avail_mem}M \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.bam \\
        --SORT_ORDER $sort_order
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
