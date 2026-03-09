process VARPIPE_TAR {

    label 'process_single'

    container 'cdc-amd/ubuntu:24.04_e2e9f67_v9'

    input:
    val(metas)
    path(files)
    val(date)

    output:
    path("*.tar.gz"), emit: varpipe_zip
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "varpipe_outputs_${date}"
    def stage_files = [
        metas.collect { "mkdir -p ${prefix}/${it.id}" },
        metas.collect { "mv ${it.filenames.join(' ')} ${prefix}/${it.id}" }
    ].flatten().join(" && ")
    """
    ${stage_files}

    tar \\
        -czh \\
        -f ${prefix}.tar.gz \\
        ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "varpipe_outputs_${date}"
    """
    touch ${prefix}.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
