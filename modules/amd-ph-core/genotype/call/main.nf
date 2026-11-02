// GENOTYPE_CALL — tabular BLAST hits in, one genotype call per query out.
//
// Deliberately knows nothing about any pathogen. What makes a call a "lineage", a "clade" or a
// "serotype" is the naming in the database it was searched against, which is an input. The only
// judgement this module makes is how much better the best subject is than the next-best DIFFERENT
// subject, and whether that margin is big enough to commit to.
//
// A typing database is dense on purpose, so the top two subjects are frequently close. Every row
// therefore carries the runner-up and the margin, and anything inside `tie_margin` is labelled
// `ambiguous` rather than silently resolved by argmax.

process GENOTYPE_CALL {
    tag "${meta.id}"
    label 'process_single'

    container 'oamd-bio-python:3.12-orb'

    input:
    tuple val(meta), path(hits)
    val genotype_opts

    output:
    tuple val(meta), path("*.genotype.tsv"), emit: calls
    // Versions go to the `versions` TOPIC rather than a versions.yml path, matching the two nf-core
    // blast modules this is used alongside. A pipeline collects every reporter with one
    // `channel.topic('versions')` instead of threading a version channel through each subworkflow.
    // The channel name on the version tuple below is required, not decorative. Without it
    // nf-core 4.1.0 lint CRASHES on this module: `_parse_output_emits` takes its unnamed-output
    // branch and appends a 3-element failure tuple to a 4-field LintResult, raising
    // `TypeError: LintResult.__init__() missing 1 required positional argument: 'file_path'`
    // and aborting the entire run. It also makes `main_nf_version_emit` balance.
    //
    // Note the parser reads EVERY non-empty line in this block, comments included, and any line
    // carrying the channel-declaration keyword is taken for a real channel. Spelling that keyword
    // out in a comment here registers a nameless output and makes `--fix` fail with
    // "list index out of range", so this note deliberately describes it instead.
    tuple val("${task.process}"), val("genotype_call"), eval("python3 --version 2>&1 | sed 's/^Python //'"), emit: versions_genotype, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Settings arrive as one Map, matching OUROBOROS_RECRUIT/POLISH, so a caller adds a knob without
    // re-ordering positional inputs. Optional with defaults: an unset key is not an error, because the
    // useful default here is "report everything and let the margin speak".
    def o          = genotype_opts ?: [:]
    def columns    = o.get('columns', 'qseqid sseqid pident length qstart qend qlen slen bitscore')
    def min_pident = o.get('min_pident', 0.0)
    def min_qcov   = o.get('min_qcov', 0.0)
    def tie_margin = o.get('tie_margin', 0.02)
    """
    genotype_call.py \\
        --hits ${hits} \\
        --sample ${meta.id} \\
        --columns '${columns}' \\
        --min-pident ${min_pident} \\
        --min-qcov ${min_qcov} \\
        --tie-margin ${tie_margin} \\
        -o ${meta.id}.genotype.tsv
    """

    stub:
    """
    printf 'sample\\tquery\\tsubject\\tpident\\tqcov\\tbitscore\\trunnerup\\trunnerup_bitscore\\tmargin\\tcall\\n' > ${meta.id}.genotype.tsv
    """
}
