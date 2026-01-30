process VARPIPE_PARSEANNOTATIONS {

    tag "${meta.id}"
    label 'process_single'

    container 'oamd-bio-python:3.12-tb_a159d0c_v2'

    input:
    tuple val(meta), path(annotation)
    path(mutationloci)

    output:
    tuple val(meta), path("*_DR_loci_Final_annotation.txt"), optional:true, emit: drLociFinalAnnotation
    tuple val(meta), path("*_full_Final_annotation.txt")   , optional:true, emit: fullFinalAnnotation
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'parse_annotation.py'

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
