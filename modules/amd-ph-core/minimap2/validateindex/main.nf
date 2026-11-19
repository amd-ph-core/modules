process MINIMAP2_VALIDATEINDEX {
    tag "$meta.id"
    label 'process_single'

    container 'cdc-amd/python:3.12-tb_716b886_v4'

    input:
    tuple val(meta), path(index), path(fasta)
    val preset

    output:
    // The index is passed through unchanged, as a link under validated/, once
    // it has been checked against the FASTA. Nextflow does not re-emit input
    // files, hence the link.
    //
    // versions.yml is a real file on the versions topic. The script is a Python
    // template, and Nextflow allows an `eval` output only with Bash scripts.
    tuple val(meta), path("validated/*"), emit: index
    path "versions.yml"                 , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'validate_minimap2_index.py'

    stub:
    """
    mkdir validated
    touch validated/${index.name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
