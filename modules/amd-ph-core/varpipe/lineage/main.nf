process VARPIPE_LINEAGE {

    tag "${meta.id}"
    label 'process_single'

    container 'cdc-amd/python:3.12-tb_716b886_v4'

    input:
    tuple val(meta), path(fullFinalAnnotation)
    path(lineages)

    output:
    // Declared as the `versions` TOPIC while still being a real file. This
    // module's script is a Python `template`, and Nextflow rejects an `eval`
    // output on any non-Bash interpreter -- "Process output of type 'eval' is
    // only allowed with Bash process scripts". Lint accepts `path "versions.yml"`
    // on a topic line for exactly this reason, so the template keeps writing the
    // file and the topic carries it.
    tuple val(meta), path("*_Lineage.txt"), emit: lineage
    tuple val(meta), path("*.lineage_report.txt"), emit: lineage_report
    path "versions.yml", emit: versions, topic: versions

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
