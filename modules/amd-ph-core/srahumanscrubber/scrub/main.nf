process SRAHUMANSCRUBBER_SCRUB {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sra-human-scrubber:2.2.1--hdfd78af_0'
        : 'biocontainers/sra-human-scrubber:2.2.1--hdfd78af_0'}"

    input:
    tuple val(meta), path(reads, stageAs: "input_reads/")
    path db

    output:
    tuple val(meta), path("*.clean*.fq.gz")         , emit: fastq
    tuple val(meta), path("*.srahumanscrubber.json"), emit: json
    tuple val(meta), path("*.srahumanscrubber.log") , emit: log
    // With -r: the identified reads as they were before masking or removal.
    tuple val(meta), path("*.removed*.fq.gz")       , emit: removed, optional: true
    // scrub.sh has no --version. The bioconda wrapper names its install directory
    // after the release (SCRUBBER_SHARE=/usr/local/share/sra-human-scrubber-2.2.1),
    // so the version is read from there. `grep .` fails on empty input, so a
    // layout without that line records "unknown" rather than an empty string.
    tuple val("${task.process}"), val('sra-human-scrubber'), eval('grep -m1 -o "sra-human-scrubber-[0-9][0-9.]*" "$(command -v scrub.sh)" | cut -d - -f 4 | grep . || echo unknown'), topic: versions, emit: versions_srahumanscrubber
    tuple val("${task.process}"), val('aligns_to'), eval('aligns_to 2>&1 | grep -o "aligns_to version [0-9.]*" | head -n 1 | cut -d " " -f 3 | grep . || echo unknown'), topic: versions, emit: versions_aligns_to

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The interface is scrub.sh's: identified reads are masked with N by default,
    // and -x removes them instead. -r also saves the identified reads, unmasked,
    // as scrub.sh -r does (cut_spots_fastq.py writes them in the same pass). No
    // other scrub.sh option is supported here.
    def unknown = args.tokenize().findAll { arg -> !(arg in ['-x', '-r']) }
    if (unknown) {
        error "SRAHUMANSCRUBBER_SCRUB: unsupported ext.args ${unknown}; only -x (remove rather than mask) and -r (save identified reads) are supported"
    }
    def save_removed = args.tokenize().contains('-r')
    def removed_flag = save_removed ? '-r removed.fastq' : ''
    def mode         = args.tokenize().contains('-x') ? 'remove' : 'mask'
    def replace_flag = mode == 'mask' ? '-n' : ''
    // Paired reads are interleaved and scrubbed as pairs (-I, scrub.sh's -s): a
    // spot is a pair, so when either mate is identified both are masked or removed,
    // and the two output files stay in step.
    def pair_flag    = meta.single_end ? '' : '-I'
    // Sorted on the full file name so mate 1 is always first (see deacon/filter).
    def sorted_reads = meta.single_end ? [reads].flatten() : reads.sort { read -> read.name }
    def in_1     = shellQuote(sorted_reads[0])
    def in_2     = meta.single_end ? '' : shellQuote(sorted_reads[1])
    def log_out  = shellQuote("${prefix}.srahumanscrubber.log")
    def json_out = shellQuote("${prefix}.srahumanscrubber.json")
    def prepare_input
    def write_output
    if (meta.single_end) {
        prepare_input = "gzip -cdf ${in_1} > input.fastq"
        write_output  = "write_reads() { compress < \$1.fastq > ${shellQuote(prefix)}.\$1.fq.gz; }"
    } else {
        prepare_input = """\
gzip -cdf ${in_1} | paste - - - - > r1.tsv
    gzip -cdf ${in_2} | paste - - - - > r2.tsv
    if [ "\$(wc -l < r1.tsv)" -ne "\$(wc -l < r2.tsv)" ]; then
        printf 'ERROR: %s and %s hold different numbers of records\\n' ${in_1} ${in_2} >&2
        exit 1
    fi
    paste r1.tsv r2.tsv | tr '\\t' '\\n' > input.fastq
    rm r1.tsv r2.tsv"""
        write_output  = """\
write_reads() {
        paste - - - - - - - - < \$1.fastq > \$1.tsv
        cut -f 1-4 \$1.tsv | tr '\\t' '\\n' | compress > ${shellQuote(prefix)}.\$1_1.fq.gz
        cut -f 5-8 \$1.tsv | tr '\\t' '\\n' | compress > ${shellQuote(prefix)}.\$1_2.fq.gz
        rm \$1.tsv
    }"""
    }
    """
    # A failure anywhere in a pipe fails the task: a gzip that dies on a truncated
    # input would otherwise be hidden by the paste after it.
    set -o pipefail

    # BGZF when the image has bgzip (as deacon/filter does), plain gzip otherwise.
    if command -v bgzip > /dev/null 2>&1; then
        compress() { bgzip --threads ${task.cpus} -c; }
    else
        compress() { gzip -n -c; }
    fi
    count_records() { echo \$(( \$(wc -l < "\$1") / 4 )); }
    # A masked record keeps its length, with every base N: count the sequence lines
    # (the second of each four) that are one or more Ns and nothing else.
    count_masked() { paste - - - - < "\$1" | cut -f 2 | { grep -c -E '^[Nn]+\$' || true; }; }

    ${prepare_input}
    reads_in=\$(count_records input.fastq)

    # scrub.sh's three steps, run directly rather than through scrub.sh: the
    # bioconda build of 2.2.1 (BioContainer 2.2.1--hdfd78af_0) rewrites scrub.sh
    # and drops the pair flag from the FASTA conversion, so with -s it numbers
    # spots per read but cuts per pair and removes the wrong reads. Passing the
    # flag to both helpers here gives the release's behaviour in any image.
    if [ "\$reads_in" -gt 0 ]; then
        fastq_to_fasta.py ${pair_flag} < input.fastq > reads.fasta
        aligns_to -db ${shellQuote(db)} -num_threads ${task.cpus} reads.fasta \\
            > identified_spots.txt 2> ${log_out}
        cut_spots_fastq.py input.fastq ${replace_flag} ${pair_flag} ${removed_flag} < identified_spots.txt \\
            > clean.fastq 2>> ${log_out}
        rm reads.fasta identified_spots.txt
    else
        # cut_spots_fastq.py writes one empty record for an empty input.
        : > clean.fastq
        echo "no reads in input" > ${log_out}
    fi

    reads_out=\$(count_records clean.fastq)
    reads_masked=\$(count_masked clean.fastq)

    ${write_output}
    write_reads clean
    if ${save_removed}; then
        # cut_spots_fastq.py creates the file only when some spot was identified.
        [ -e removed.fastq ] || : > removed.fastq
        write_reads removed
    fi
    rm -f input.fastq clean.fastq removed.fastq

    # Counts are records, not pairs, as in the deacon and hostile reports, so a
    # paired-end sample's numbers are twice its pair counts.
    printf '{\\n  "tool": "sra-human-scrubber",\\n  "database": %s,\\n  "mode": "%s",\\n  "paired": %s,\\n  "reads_in": %d,\\n  "reads_out": %d,\\n  "reads_removed": %d,\\n  "reads_masked": %d\\n}\\n' \\
        ${shellQuote(groovy.json.JsonOutput.toJson(db.name.toString()))} "${mode}" "${!meta.single_end}" \\
        "\$reads_in" "\$reads_out" "\$(( reads_in - reads_out ))" "\$reads_masked" \\
        > ${json_out}
    """

    stub:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def mode    = args.tokenize().contains('-x') ? 'remove' : 'mask'
    def kinds   = args.tokenize().contains('-r') ? ['clean', 'removed'] : ['clean']
    def outputs = kinds.collectMany { kind -> meta.single_end ? ["${prefix}.${kind}.fq.gz"] : ["${prefix}.${kind}_1.fq.gz", "${prefix}.${kind}_2.fq.gz"] }
    """
    for f in ${outputs.collect { o -> shellQuote(o) }.join(' ')}; do
        echo "" | gzip > "\$f"
    done
    echo "srahumanscrubber stub" > ${shellQuote("${prefix}.srahumanscrubber.log")}
    printf '{\\n  "tool": "sra-human-scrubber",\\n  "database": %s,\\n  "mode": "%s",\\n  "paired": %s,\\n  "reads_in": 0,\\n  "reads_out": 0,\\n  "reads_removed": 0,\\n  "reads_masked": 0\\n}\\n' \\
        ${shellQuote(groovy.json.JsonOutput.toJson(db.name.toString()))} "${mode}" "${!meta.single_end}" > ${shellQuote("${prefix}.srahumanscrubber.json")}
    """
}

//
// Single-quote a path for the shell, closing and reopening around any embedded
// quote, so a staged name or prefix containing whitespace stays one argument.
//
def shellQuote(path) {
    return "'" + path.toString().replace("'", "'\\''") + "'"
}
