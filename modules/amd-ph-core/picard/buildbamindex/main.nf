process PICARD_BUILDBAMINDEX {

    tag "${meta.id}"
    label 'process_medium'

    container 'cdc-amd/picard:3.4.0_ac49616_v8'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bai"), emit: bai
    tuple val("${task.process}"), val('picard'), eval('picard CreateSequenceDictionary --version 2>&1 | grep -o "Version:.*" | cut -f2- -d:'), emit: versions_picard, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
        picard BuildBamIndex \\
            INPUT=${input} \\
            VALIDATION_STRINGENCY=${args}
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bai
    """
}
