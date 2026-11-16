process SAMTOOLS_ALIGNMENT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'cdc-amd/samtools:1.22.1_cc13003_v8'

    input:
    // The reference travels in the same tuple as the alignment rather than on a
    // channel of its own. A pipeline that aligns every sample against several
    // references produces one BAM per sample/reference pair, and a separate
    // reference channel would be paired with the BAMs by arrival order -- a BAM
    // could be handed another reference's FASTA. Pass [] for fasta and fai when
    // there is no reference; it is only required for CRAM input.
    tuple val(meta), path(input), path(index), path(fasta), path(fai)
    val run_stats

    output:
    tuple val(meta), path("*.coverage.txt")   , emit: coverage
    tuple val(meta), path("*.depth.tsv")      , emit: depth
    tuple val(meta), path("*.stats")          , emit: stats    , optional: true
    tuple val(meta), path("*.flagstat")       , emit: flagstat , optional: true
    tuple val(meta), path("*.idxstats")       , emit: idxstats , optional: true
    tuple val(meta), path("*.samtools_qc.tar"), emit: qc_bundle, optional: true
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | sed -n "1s/^samtools //p"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Every samtools step reads the same staged alignment inside one task, so
    // the alignment is staged once instead of once per step. nf-core lint
    // allows only args/args2/args3 as ext keys, so they go to the three steps
    // with options worth tuning; flagstat and idxstats take none.
    def args          = task.ext.args  ?: ''   // samtools stats
    def args2         = task.ext.args2 ?: ''   // samtools coverage
    def args3         = task.ext.args3 ?: ''   // samtools depth
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def reference     = fasta ? "--reference ${fasta}" : ''
    def extra_threads = Math.max(0, task.cpus - 1)
    // The QC files are also packed into one archive so a downstream report
    // stages one file per alignment instead of three. The archive keeps the
    // layout MultiQC expects: flagstat/ and samtools/.
    def stats_command = run_stats ? """
    samtools \\
        stats \\
        --threads ${task.cpus} \\
        ${args} \\
        ${reference} \\
        ${input} \\
        > ${prefix}.stats

    samtools \\
        flagstat \\
        --threads ${task.cpus} \\
        ${input} \\
        > ${prefix}.flagstat

    samtools \\
        idxstats \\
        --threads ${extra_threads} \\
        ${input} \\
        > ${prefix}.idxstats

    mkdir -p qc_bundle/flagstat qc_bundle/samtools
    cp ${prefix}.flagstat qc_bundle/flagstat/
    cp ${prefix}.stats ${prefix}.idxstats qc_bundle/samtools/
    tar \\
        -c \\
        -f ${prefix}.samtools_qc.tar \\
        -C qc_bundle .
    """ : ''
    """
    ${stats_command}

    samtools \\
        coverage \\
        ${args2} \\
        ${reference} \\
        -o ${prefix}.coverage.txt \\
        ${input}

    samtools \\
        depth \\
        --threads ${extra_threads} \\
        ${args3} \\
        -o ${prefix}.depth.tsv \\
        ${input}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stats_command = run_stats ? """
    touch ${prefix}.stats
    touch ${prefix}.flagstat
    touch ${prefix}.idxstats
    mkdir -p qc_bundle/flagstat qc_bundle/samtools
    cp ${prefix}.flagstat qc_bundle/flagstat/
    cp ${prefix}.stats ${prefix}.idxstats qc_bundle/samtools/
    tar -c -f ${prefix}.samtools_qc.tar -C qc_bundle .
    """ : ''
    """
    ${stats_command}
    touch ${prefix}.coverage.txt
    touch ${prefix}.depth.tsv
    """
}
