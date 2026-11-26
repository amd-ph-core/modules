process FALCO {
    tag "${meta.id}"
    // process_low rather than process_single: paired-end input runs one falco
    // per mate concurrently, so the task wants two cores. Single-end input
    // leaves the second idle, which is the cheaper mistake — the alternative is
    // a one-core reservation that paired-end input immediately doubles.
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/falco:2.0.3--h3be2455_0'
        : 'biocontainers/falco:2.0.3--h3be2455_0'}"

    input:
    tuple val(meta), path(reads)

    // Optional: a mate with no reads gets no report (falco 2.x fails on an empty
    // file), so an empty sample emits nothing here rather than failing the task.
    output:
    tuple val(meta), path("*_fastqc_report.html"), emit: html   , optional: true
    tuple val(meta), path("*_fastqc_data.txt")   , emit: txt    , optional: true
    tuple val(meta), path("*_summary.txt")       , emit: summary, optional: true
    tuple val("${task.process}"), val('falco'), eval("falco --version | head -n 1"), topic: versions, emit: versions_falco

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Each mate is named from its position in the staged list, `<prefix>_R1` and
    // `<prefix>_R2` (or just `<prefix>` single-end), not from its filename, so the
    // same sample reports under the same name whatever the sequencer or an
    // upstream step called the file.
    def files   = reads.toList()
    def links   = []
    def runs    = []
    def renames = []
    def outdirs = []
    files.eachWithIndex { read, i ->
        def name = files.size() == 1 ? prefix : "${prefix}_R${i + 1}"
        def link = "${name}${fastqExtension(read)}"
        def out  = "falco_${name}"
        links   << "ln -s \"\$(readlink -f ${shellQuote(read)})\" ${shellQuote(link)}"
        runs    << """if has_reads ${shellQuote(link)}; then
        falco ${args} --threads 1 --output ${shellQuote(out)} ${shellQuote(link)} &
        pids="\$pids \$!"
    else
        echo FALCO: ${shellQuote(name)}: no reads, no report >&2
    fi"""
        renames << """if [ -d ${shellQuote("${out}/${name}")} ]; then
        mv ${shellQuote("${out}/${name}/fastqc_data.txt")} ${shellQuote("${name}_fastqc_data.txt")}
        mv ${shellQuote("${out}/${name}/summary.txt")} ${shellQuote("${name}_summary.txt")}
        mv ${shellQuote("${out}/${name}/fastqc_report.html")} ${shellQuote("${name}_fastqc_report.html")}
    fi"""
        outdirs << shellQuote(out)
    }
    // falco is handed a LINK named for the mate rather than the staged file. falco
    // records the name of the file it was given as `Filename`, and MultiQC's FastQC
    // parser takes the sample name from that embedded value, not from the report's
    // path; the extension is kept so falco still detects compression.
    //
    // Each mate runs as its own single-threaded falco, concurrently. falco 2.x can
    // take both mates in one call on --threads, but then its per-tile quality
    // figures depend on how the work is split between threads: measured on 2.0.3,
    // two runs of the same pair at --threads 2 gave different fastqc_data.txt, and
    // one in six differed under CPU contention, while --threads 1 was identical
    // every time. A QC report that changes from run to run on the same reads is not
    // one to publish, so the parallelism stays per mate.
    //
    // falco 2.x writes OUTDIR/<input name>/{fastqc_data.txt,summary.txt,fastqc_report.html};
    // the outputs are renamed to the names MultiQC's FastQC module matches
    // (`*fastqc_data.txt`). It reads its adapter, contaminant and limits files from
    // its install prefix's share/falco.
    //
    // Each PID is waited on individually so its status is checked: a bare `wait`
    // returns zero however the children exited, which would leave a missing or
    // truncated report behind a task that reported success.
    //
    // A mate with no reads is skipped: falco 2.x fails on an empty file (2.0.4: exit 1,
    // "failed to read header", gzipped; a floating-point exception, plain), which would
    // stop a whole run over one empty sample. has_reads decompresses (or passes through)
    // only the first byte; `|| true` keeps head's early close from failing the pipeline.
    // A 0-byte file is empty (some gzip builds reject it as "invalid magic" even with -f).
    // Any other file that yields no bytes is decompressed in full, with errors kept, so a
    // damaged one (a gzip truncated to its header also yields nothing) fails the task
    // instead of passing as empty. `exit` is explicit: errexit does not apply inside an
    // `if` condition.
    """
    has_reads() {
        if [ ! -s "\$1" ]; then return 1; fi
        if [ -n "\$(gzip -cdf -- "\$1" 2>/dev/null | head -c 1 || true)" ]; then return 0; fi
        gzip -cdf -- "\$1" > /dev/null || { echo "FALCO: \$1: not a readable FASTQ" >&2; exit 1; }
        return 1
    }

    ${links.join('\n    ')}

    pids=''
    ${runs.join('\n    ')}
    for pid in \$pids; do
        wait "\$pid"
    done

    ${renames.join('\n    ')}
    rm -rf ${outdirs.join(' ')}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (reads.toList().size() == 1) {
        """
        touch ${shellQuote("${prefix}_fastqc_data.txt")}
        touch ${shellQuote("${prefix}_summary.txt")}
        touch ${shellQuote("${prefix}_fastqc_report.html")}
        """
    }
    else {
        """
        for mate in 1 2; do
            read_base=${shellQuote(prefix)}"_R\${mate}"
            touch "\${read_base}_fastqc_data.txt"
            touch "\${read_base}_summary.txt"
            touch "\${read_base}_fastqc_report.html"
        done
        """
    }
}

//
// The extension falco needs to see to detect compression: falco 2.x reads plain or
// gzip/BGZF FASTQ and decides by the name.
//
def fastqExtension(read) {
    return read.name.endsWith('.gz') ? '.fastq.gz' : '.fastq'
}

//
// Single-quote a path for the shell, closing and reopening around any embedded
// quote, so a file name or prefix with an apostrophe or a space stays one argument.
//
def shellQuote(path) {
    return "'" + path.toString().replace("'", "'\\''") + "'"
}
