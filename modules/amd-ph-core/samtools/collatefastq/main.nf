process SAMTOOLS_COLLATEFASTQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.25--h9dcdb79_0'
        : 'biocontainers/samtools:1.25--h9dcdb79_0'}"

    input:
    // Staged into a subdirectory so the output globs below can never match an input.
    tuple val(meta), path(input, stageAs: "input/*")
    // Reference, for CRAM only. Pass [ [:], [], [] ] for BAM or SAM.
    tuple val(meta2), path(fasta), path(fai)

    output:
    // All four files are always written, holding no reads where there are none, so a
    // caller never has to test for a file's presence. Which of them hold reads is in
    // `read_counts`, as values: that is how a caller learns the input's layout.
    tuple val(meta), path("*_{1,2}.fq.gz"),   emit: fastq
    tuple val(meta), path("*_other.fq.gz"),     emit: fastq_other
    tuple val(meta), path("*_singleton.fq.gz"), emit: fastq_singleton
    // Reads in each file, in the order read 1, read 2, other, singleton.
    //   paired-end input:  read_1 = read_2 > 0; other + singleton are orphans
    //   single-end input:  read_1 = read_2 = 0; the reads are in other
    tuple val(meta), eval("cat read_counts/1"), eval("cat read_counts/2"), eval("cat read_counts/other"), eval("cat read_counts/singleton"), emit: read_counts
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${shellQuote(fasta)}" : ""
    // collate's fast mode outputs only records that are paired: on a single-end BAM it
    // writes nothing and exits 0. Measured on samtools 1.22.1 -- 0 of 2,000 reads.
    // -f is refused on its own and clustered with collate's other no-value flags
    // (-O, -u), as in `-uf`, which samtools reads as `-u -f`.
    if (args.tokenize().any { token -> token ==~ /-[Ouf]*f[Ouf]*/ }) {
        error "SAMTOOLS_COLLATEFASTQ: `samtools collate -f` drops every unpaired read, so a single-end input comes out empty. Remove -f from ext.args."
    }
    """
    # A failure early in a pipe must fail the task: a collate that dies after the
    # header would otherwise leave fastq writing a truncated conversion that exits 0.
    set -o pipefail

    # Every input is collated, so mates are adjacent whatever the sort order: a
    # coordinate-sorted BAM fed straight to `samtools fastq` sends the mates that
    # are not adjacent to the singleton file. `samtools fastq` drops secondary and
    # supplementary records by default (-F 0x900); set ext.args2 to change that.
    samtools collate \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reference} \\
        -O \\
        ${shellQuote(input)} \\
        collate_tmp |
    samtools fastq \\
        ${args2} \\
        --threads ${task.cpus} \\
        ${reference} \\
        -1 ${shellQuote("${prefix}_1.fq.gz")} \\
        -2 ${shellQuote("${prefix}_2.fq.gz")} \\
        -0 ${shellQuote("${prefix}_other.fq.gz")} \\
        -s ${shellQuote("${prefix}_singleton.fq.gz")}

    # Count the reads written to each file. `samtools view -c` cannot do it: it
    # fails on a FASTQ that holds no reads instead of returning 0.
    mkdir read_counts
    for part in 1 2 other singleton; do
        lines=\$(gzip -cd ${shellQuote(prefix)}"_\${part}.fq.gz" | wc -l)
        echo \$(( lines / 4 )) > "read_counts/\${part}"
    done
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir read_counts
    for part in 1 2 other singleton; do
        echo | gzip > "${prefix}_\${part}.fq.gz"
        echo 0 > "read_counts/\${part}"
    done
    """
}

//
// Single-quote a path for the shell, closing and reopening around any embedded
// quote, so a staged name or prefix containing whitespace stays one argument.
//
def shellQuote(path) {
    return "'" + path.toString().replace("'", "'\\''") + "'"
}
