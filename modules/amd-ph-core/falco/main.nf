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

    output:
    tuple val(meta), path("*_fastqc_report.html"), emit: html
    tuple val(meta), path("*_fastqc_data.txt")   , emit: txt
    tuple val(meta), path("*_summary.txt")       , emit: summary
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
        links   << "ln -s \"\$(readlink -f '${read}')\" '${link}'"
        runs    << "falco ${args} --threads 1 --output 'falco_${name}' '${link}' &\n    pids=\"\$pids \$!\""
        renames << "mv 'falco_${name}/${name}/fastqc_data.txt'    '${name}_fastqc_data.txt'"
        renames << "mv 'falco_${name}/${name}/summary.txt'        '${name}_summary.txt'"
        renames << "mv 'falco_${name}/${name}/fastqc_report.html' '${name}_fastqc_report.html'"
        outdirs << "'falco_${name}'"
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
    """
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
        touch ${prefix}_fastqc_data.txt
        touch ${prefix}_summary.txt
        touch ${prefix}_fastqc_report.html
        """
    }
    else {
        """
        for mate in 1 2; do
            read_base="${prefix}_R\${mate}"
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
