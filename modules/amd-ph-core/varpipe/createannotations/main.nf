process VARPIPE_CREATEANNOTATIONS {

    tag "${meta.id}"
    label 'process_single'

    container 'cdc-amd/python:3.12-tb_716b886_v4'

    input:
    tuple val(meta), path(annotation)

    output:
    // Declared as the `versions` TOPIC while still being a real file. This
    // module's script is a Python `template`, and Nextflow rejects an `eval`
    // output on any non-Bash interpreter -- "Process output of type 'eval' is
    // only allowed with Bash process scripts". Lint accepts `path "versions.yml"`
    // on a topic line for exactly this reason, so the template keeps writing the
    // file and the topic carries it.
    tuple val(meta), path("*_DR_loci_annotation.txt")      , optional:true, emit: drLociAnnotation
    tuple val(meta), path("*_full_annotation.txt")         , optional:true, emit: fullAnnotation
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'create_annotation.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_DR_loci_Final_annotation.txt
    touch ${prefix}_full_Final_annotation.txt
    touch ${prefix}_DR_loci_annotation.txt
    touch ${prefix}_full_annotation.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
