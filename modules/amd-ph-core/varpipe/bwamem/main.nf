process VARPIPE_BWAMEM {
    tag "$meta.id"
    label 'process_high'

    container 'cdc-amd/bwa:0.7.19_639f31d_v9'

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index), path(ref_dict), path(ref_fai), path(fasta)

    output:
    tuple val(meta), path("*.sam")  , emit: sam
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rg_line = "-R '@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:ILLUMINA'"
    def args = task.ext.args ? "${task.ext.args} ${rg_line}" : "${rg_line}"
    """
    cp bwa/* .
    INDEX=`find -L bwa/ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa mem \\
        $args \\
        \$INDEX \\
        ${reads} \\
        > ${prefix}.sam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
