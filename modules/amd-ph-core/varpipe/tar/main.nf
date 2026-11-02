process VARPIPE_TAR {

    label 'process_single'

    container 'cdc-amd/ubuntu:24.04_e2e9f67_v9'

    input:
    val(metas)
    path(files)
    val(date)

    output:
    path("*.tar.gz"), emit: varpipe_zip
    tuple val("${task.process}"), val('tar'), eval('tar --version | head -1 | sed "s/^tar (GNU tar) //"'), emit: versions_tar, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "varpipe_outputs_${date}"
    // `echo "" | gzip >` rather than `touch`. A touched `.tar.gz` is a zero-byte
    // file, not a valid gzip stream, so anything that opens it as one -- nf-test's
    // gzip parser included -- fails on a file that merely looks present. The
    // linter accepts this exact form only: its regex is
    // `echo\s+""\s*\|\s*gzip\s*>\s*.*\.gz$`, so `gzip -n` is rejected even though
    // it is the more reproducible flag.
    """
    echo "" | gzip > ${prefix}.tar.gz
    """
}
