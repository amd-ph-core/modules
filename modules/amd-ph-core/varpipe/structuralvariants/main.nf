process VARPIPE_STRUCTURALVARIANTS {

    tag "${meta.id}"
    label 'process_single'

    container 'cdc-amd/python:3.12-tb_716b886_v4'

    input:
    tuple val(meta), path(qc_status), path(coverage), path(unmapped), path(mapped)
    tuple val(meta2), path(genome_region_coverage)
    path(bedstruct)

    output:
    tuple val(meta), path("*_structural_variants.txt"), emit: structural_variants, optional: true
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'structvar_parser.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_stats.txt
    touch ${prefix}_target_region_coverage.txt
    touch ${prefix}_genome_region_coverage.txt
    touch ${prefix}_structural_variants.txt
    touch .qc_passed.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
