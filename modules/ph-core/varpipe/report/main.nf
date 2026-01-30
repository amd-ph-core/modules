process VARPIPE_REPORT {

    tag "${meta.id}"
    label 'process_single'

    container 'oamd-bio-python:3.12-tb_a159d0c_v2'

    input:
    tuple val(meta), path(statsFile), path(targetRegionCoverage), path(drLociFinalAnnotation)

    output:
    tuple val(meta), path("summary.txt"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'create_report.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
