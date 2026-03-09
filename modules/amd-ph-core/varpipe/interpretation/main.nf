process VARPIPE_INTERPRETATION {

    tag "${meta.id}"
    label 'process_single'

    container 'cdc-amd/python:3.12-tb_716b886_v4'

    input:
    tuple val(meta), path(summary), path(structuralVariants), path(drLociAnnotation), path(targetRegionCoverage)
    path(reported)

    output:
    tuple val(meta), path("*_summary.txt"), emit: summary
    tuple val(meta), path("*_interpretation.txt"), emit: interpretation
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'interpret.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_summary.txt
    touch ${prefix}_interpretation.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
