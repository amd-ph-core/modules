process PICARD_CREATESEQUENCEDICTIONARY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'oamd-bio-picard:3.4.0_3935237_v4'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.dict"), emit: reference_dict
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CreateSequenceDictionary] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        CreateSequenceDictionary  \\
        $args \\
        --REFERENCE $fasta \\
        --OUTPUT ${prefix}.dict

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CreateSequenceDictionary --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dict

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CreateSequenceDictionary --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

}
