process PICARD_SAMFORMATCONVERTER {

    tag "${meta.id}"
    label 'process_medium'

    container 'cdc-amd/picard:3.4.0_ac49616_v8'

    input:
    tuple val(meta), path(alnSam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_suffix = task.ext.output_suffix ?: 'bam'
    """
    picard SamFormatConverter \\
        INPUT=${alnSam} \\
        VALIDATION_STRINGENCY=${args} \\
        OUTPUT=${prefix}.${output_suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard SamFormatConverter --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard SamFormatConverter --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
