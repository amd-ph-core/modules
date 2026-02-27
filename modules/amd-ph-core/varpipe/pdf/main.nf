process VARPIPE_PDF {

    tag "${meta.id}"
    label 'process_single'

    container 'oamd-bio-python:3.12-tb_a159d0c_v2'

    input:
    tuple val(meta), path(summary)

    output:
    tuple val(meta), path("*_report.pdf"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'pdf_print.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
