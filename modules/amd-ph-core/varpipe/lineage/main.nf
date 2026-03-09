process VARPIPE_LINEAGE {

    tag "${meta.id}"
    label 'process_single'

    container 'cdc-amd/python:3.12-tb_716b886_v4'

    input:
    tuple val(meta), path(fullFinalAnnotation)
    path(lineages)

    output:
    tuple val(meta), path("*_Lineage.txt"), emit: lineage
    tuple val(meta), path("*.lineage_report.txt"), emit: lineage_report
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'lineage_parser.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_Lineage.txt
    touch ${prefix}.lineage_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
