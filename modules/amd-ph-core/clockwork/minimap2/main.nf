process CLOCKWORK_MINIMAP2 {
    tag "$meta.id"
    label 'process_medium'

    container 'cdc-amd/minimap2:2.30_3457b81_v6'

    input:
    tuple val(meta), path(reads)
    path(reference)

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path("versions.yml")          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def split_prefix = "--split-prefix ${prefix}"
    def read_group = "-R '@RG\\tLB:LIB\\tID:1\\tSM:${prefix}'"
    def contam_sam = "${prefix}.contam.sam"
    def sam_filter_output = "-a | samtools view -h -e '!(flag.secondary || flag.supplementary)' -o ${contam_sam}"

    """
    minimap2 \\
        -x sr \\
        -t $task.cpus \\
        $read_group \\
        $split_prefix \\
        $reference \\
        $reads \\
        $sam_filter_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def contam_sam = "${prefix}.contam.sam"

    """
    touch $contam_sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
