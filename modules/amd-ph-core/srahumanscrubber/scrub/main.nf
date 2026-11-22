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
    // scrub.sh has no --version. The bioconda wrapper names its install directory
    // after the release (SCRUBBER_SHARE=/usr/local/share/sra-human-scrubber-2.2.1),
    // so the version is read from there. `grep .` fails on empty input, so a
    // layout without that line records "unknown" rather than an empty string.
    tuple val("${task.process}"), val('sra-human-scrubber'), eval('grep -m1 -o "sra-human-scrubber-[0-9][0-9.]*" "$(command -v scrub.sh)" | sed "s/^sra-human-scrubber-//" | grep . || echo unknown'), topic: versions, emit: versions_srahumanscrubber
    tuple val("${task.process}"), val('aligns_to'), eval('aligns_to 2>&1 | sed -n "s/.*aligns_to version //p" | head -n 1 | grep . || echo unknown'), topic: versions, emit: versions_aligns_to

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // scrub.sh masks identified reads with N by default; -x removes them instead.
    def mode   = (args =~ /(^|\s)-x(\s|$)/) ? 'remove' : 'mask'
    // Sorted on the full file name so mate 1 is always first (see deacon/filter).
    def sorted_reads = meta.single_end ? [reads].flatten() : reads.sort { read -> read.name }
    def prepare_input
    def write_output
    if (meta.single_end) {
        prepare_input = "gzip -cdf ${sorted_reads[0]} > input.fastq"
        write_output  = "compress < clean.fastq > ${prefix}.clean.fq.gz"
    } else {
        // scrub.sh reads one file. The mates are interleaved so that -s treats
        // each pair as a unit: when either mate is identified as human, both are
        // masked or removed, and the output stays pairwise in step.
        prepare_input = """\
gzip -cdf ${sorted_reads[0]} | paste - - - - > r1.tsv
    gzip -cdf ${sorted_reads[1]} | paste - - - - > r2.tsv
    if [ "\$(wc -l < r1.tsv)" -ne "\$(wc -l < r2.tsv)" ]; then
        echo "ERROR: ${sorted_reads[0]} and ${sorted_reads[1]} hold different numbers of records" >&2
        exit 1
    fi
    paste r1.tsv r2.tsv | tr '\\t' '\\n' > input.fastq
    rm r1.tsv r2.tsv"""
        write_output  = """\
paste - - - - - - - - < clean.fastq > clean.tsv
    cut -f 1-4 clean.tsv | tr '\\t' '\\n' | compress > ${prefix}.clean_1.fq.gz
    cut -f 5-8 clean.tsv | tr '\\t' '\\n' | compress > ${prefix}.clean_2.fq.gz
    rm clean.tsv"""
    }
    def pair_flag = meta.single_end ? '' : '-s'
    """
    # scrub.sh writes a FASTA copy of the whole input to a mktemp directory, so
    # keep it in the task directory rather than a possibly small /tmp.
    export TMPDIR=\$PWD/tmp
    mkdir -p \$TMPDIR

    # BGZF when the image has bgzip (as deacon/filter does), plain gzip otherwise.
    if command -v bgzip > /dev/null 2>&1; then
        compress() { bgzip --threads ${task.cpus} -c; }
    else
        compress() { gzip -n -c; }
    fi
    count_records() { awk 'END { print int(NR / 4) }'; }
    # A masked record keeps its length, with every base N.
    count_masked() { awk 'NR % 4 == 2 && length(\$0) > 0 && \$0 !~ /[^Nn]/ { n++ } END { print n + 0 }'; }

    ${prepare_input}
    reads_in=\$(count_records < input.fastq)

    scrub.sh \\
        ${args} \\
        ${pair_flag} \\
        -p ${task.cpus} \\
        -d ${db} \\
        -i input.fastq \\
        -o clean.fastq \\
        2> ${prefix}.srahumanscrubber.log \\
        || { cat ${prefix}.srahumanscrubber.log >&2; exit 1; }

    reads_out=\$(count_records < clean.fastq)
    reads_masked=\$(count_masked < clean.fastq)

    ${write_output}
    rm -rf input.fastq clean.fastq \$TMPDIR

    # Counts are records, not pairs, as in the deacon and hostile reports, so a
    # paired-end sample's numbers are twice its pair counts.
    printf '{\\n  "tool": "sra-human-scrubber",\\n  "database": "%s",\\n  "mode": "%s",\\n  "paired": %s,\\n  "reads_in": %d,\\n  "reads_out": %d,\\n  "reads_removed": %d,\\n  "reads_masked": %d\\n}\\n' \\
        "${db.name}" "${mode}" "${!meta.single_end}" \\
        "\$reads_in" "\$reads_out" "\$(( reads_in - reads_out ))" "\$reads_masked" \\
        > ${prefix}.srahumanscrubber.json
    """

    stub:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def mode    = (args =~ /(^|\s)-x(\s|$)/) ? 'remove' : 'mask'
    def outputs = meta.single_end ? ["${prefix}.clean.fq.gz"] : ["${prefix}.clean_1.fq.gz", "${prefix}.clean_2.fq.gz"]
    """
    for f in ${outputs.join(' ')}; do
        echo "" | gzip > \$f
    done
    echo "srahumanscrubber stub" > ${prefix}.srahumanscrubber.log
    printf '{\\n  "tool": "sra-human-scrubber",\\n  "database": "%s",\\n  "mode": "%s",\\n  "paired": %s,\\n  "reads_in": 0,\\n  "reads_out": 0,\\n  "reads_removed": 0,\\n  "reads_masked": 0\\n}\\n' \\
        "${db.name}" "${mode}" "${!meta.single_end}" > ${prefix}.srahumanscrubber.json
    """
}
