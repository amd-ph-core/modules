process IRMA {
    tag "$meta.id"
    label 'process_high'

    container 'oamd-bio-cdcgov-irma:v1.3.0-rc1_93ea5ec_v1'

    input:
    tuple val(meta), path(reads)
    val irma_module

    output:
    tuple val(meta), path ("*/amended_consensus/*.fa"), emit: amended_consensus, optional:true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    IRMA \\
        ${irma_module} \\
        ${reads} \\
        ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma: \$(echo \$(IRMA | grep -o "v[0-9][^ ]*" | cut -c 2-))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/amended_consensus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma: \$(echo \$(IRMA | grep -o "v[0-9][^ ]*" | cut -c 2-))
    END_VERSIONS
    """
}
