process VARPIPE_STRUCTURALVARIANTS {

    tag "${meta.id}"
    label 'process_single'

    container 'cdc-amd/python:3.12-tb_716b886_v4'

    input:
    tuple val(meta), path(qc_status), path(coverage), path(unmapped), path(mapped)
    tuple val(meta2), path(genome_region_coverage)
    path(bedstruct)

    output:
    // Declared as the `versions` TOPIC while still being a real file. This
    // module's script is a Python `template`, and Nextflow rejects an `eval`
    // output on any non-Bash interpreter -- "Process output of type 'eval' is
    // only allowed with Bash process scripts". Lint accepts `path "versions.yml"`
    // on a topic line for exactly this reason, so the template keeps writing the
    // file and the topic carries it.
    tuple val(meta), path("*_structural_variants.txt"), emit: structural_variants, optional: true
    path "versions.yml", emit: versions, topic: versions

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
