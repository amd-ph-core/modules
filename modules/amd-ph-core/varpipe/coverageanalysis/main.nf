process VARPIPE_COVERAGEANALYSIS {

    tag "${meta.id}"
    label 'process_single'

    container 'cdc-amd/python:3.12-tb_716b886_v4'

    input:
    tuple val(meta), path(coverage), path(unmapped), path(mapped)
    path(amp) // gene regions (single BED file)
    path(features) // feature regions (multiple BED files)
    val(pipeline_version)

    output:
    // Declared as the `versions` TOPIC while still being a real file. This
    // module's script is a Python `template`, and Nextflow rejects an `eval`
    // output on any non-Bash interpreter -- "Process output of type 'eval' is
    // only allowed with Bash process scripts". Lint accepts `path "versions.yml"`
    // on a topic line for exactly this reason, so the template keeps writing the
    // file and the topic carries it.
    tuple val(meta), path("*_stats.txt"), emit: stats
    tuple val(meta), path("*_target_region_coverage.txt"), emit: target_region_coverage, optional: true
    tuple val(meta), path("*_genome_region_coverage.txt"), emit: genome_region_coverage, optional: true
    tuple val(meta), path("*_qc_log.txt"), emit: qc_failed, optional: true
    tuple val(meta), path(".qc_passed.txt"), emit: qc_passed, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'coverage_stats.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_stats.txt
    touch ${prefix}_target_region_coverage.txt
    touch ${prefix}_genome_region_coverage.txt
    touch .qc_passed.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
