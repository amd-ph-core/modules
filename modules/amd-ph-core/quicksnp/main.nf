process QUICKSNP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "us-cdcgov/cdc-amd/quicksnp:1.0.1_1709e03_v2"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.nwk"), emit: quicksnp_tree
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    QuickSNP.py \\
        --dm ${tsv} \\
        --outtree quicksnp_tree.nwk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch quicksnp_tree.nwk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
