process SRAHUMANSCRUBBER_BUILDDB {
    tag "${meta.id}"
    label 'process_high'

    // The STAT database tools (build_index_of_each_file, merge_db, subtract_db,
    // aligns_to) come from ncbi/ngs-tools (tools/tax) and have no conda package or
    // BioContainer: the bioconda sra-human-scrubber package ships aligns_to only.
    // A private sra-human-scrubber image carries them alongside the scrubber, built
    // from the same aligns_to, so one image serves this module and scrub.
    // Consumers map the name, and set the registry and the exact tag, in their own
    // config.
    container 'oamd-bio-sra-human-scrubber:2.2.1_a0267a0_v3'

    input:
    tuple val(meta), path(fasta, stageAs: "host/*")
    path(protect, stageAs: "protect/*")
    val window

    output:
    tuple val(meta), path("*.db")                    , emit: db
    tuple val(meta), path("*.builddb.json")          , emit: json
    tuple val("${task.process}"), val('aligns_to'), eval('aligns_to 2>&1 | grep -o "aligns_to version [0-9.]*" | head -n 1 | cut -d " " -f 3 | grep . || echo unknown'), topic: versions, emit: versions_aligns_to
    tuple val("${task.process}"), val('build_index_of_each_file'), eval('build_index_of_each_file 2>&1 | grep -o "build_index_of_each_file version [0-9.]*" | head -n 1 | cut -d " " -f 3 | grep . || echo unknown'), topic: versions, emit: versions_build_index_of_each_file

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def protect_set = protect ? [protect].flatten() : []
    def kmer_len    = 32
    def db_out      = shellQuote("${prefix}.db")
    def json_out    = shellQuote("${prefix}.builddb.json")
    def protect_in  = protect_set.collect { f -> shellQuote(f) }.join(' ')
    """
    # A failure anywhere in a pipe fails the task: a truncated compressed FASTA
    # would otherwise build an incomplete database without a word.
    set -o pipefail

    # A STAT .db is a header (format version, k) and a sorted set of canonical
    # k-mers. From each consecutive window of a sequence, build_index_of_each_file
    # keeps the k-mer with the smallest hash: window 1 keeps every k-mer, window 32
    # one in 32 (NCBI's human_filter.db is consistent with 32). aligns_to tests
    # every k-mer of a read, so an error-free host read of at least 2 x window +
    # k - 1 bases always contains a kept k-mer.
    db_kmers() { echo \$(( (\$(stat -c %s "\$1") - 24) / 8 )); }

    # Every intermediate lives in work/, and only the final database is written
    # beside it, so no prefix (`host`, `shared`, ...) can collide with one.
    mkdir -p work/plain

    # One plain, uppercase FASTA per input: the tools do not read gzip, and they
    # silently skip lowercase bases, so a soft-masked genome (repeats in
    # lowercase, about half of T2T-CHM13) would otherwise lose those k-mers.
    # Headers are uppercased too, which is harmless: the database holds k-mers only.
    upper() { gzip -cdf "\$1" | tr '[:lower:]' '[:upper:]'; }
    # Named after the full staged name, so ref.fa and ref.fa.gz stay two files.
    for f in host/*; do
        upper "\$f" > "work/plain/\$(basename "\$f").plain.fa"
    done
    : > work/host.list
    host_bytes=0
    for f in work/plain/*; do
        size=\$(stat -c %s "\$f")
        host_bytes=\$(( host_bytes + size ))
        printf '%s\\t%s\\n' "\$size" "\$PWD/\$f" >> work/host.list
    done
    build_index_of_each_file work/host.list ${window} ${kmer_len} .w${window}.db
    # merge_db sizes its hash table up front: about two slots per expected k-mer.
    expected=\$(( 2 * host_bytes / ${window} + 1000 ))
    merge_db work/host.list .w${window}.db work/host.db "\$expected"
    host_kmers=\$(db_kmers work/host.db)

    shared_kmers=0
    if ${protect_set ? 'true' : 'false'}; then
        # Protecting sequences: every host k-mer found in them is dropped, so a read
        # from them cannot be identified as host through a shared k-mer. aligns_to
        # streams the protected FASTAs against the host set and prints each host
        # k-mer it meets; that list is rebuilt into a set at window 1 (which also
        # puts each k-mer in canonical orientation) and subtracted.
        : > work/shared.txt
        for f in ${protect_in}; do
            upper "\$f" > work/protect.fa
            aligns_to -db work/host.db -num_threads ${task.cpus} -print_kmers_only work/protect.fa >> work/shared.txt
        done
        # One FASTA record per distinct k-mer, named k1, k2, ...
        sort -u work/shared.txt > work/shared.sorted
        n_shared=\$(wc -l < work/shared.sorted)
        if [ "\$n_shared" -gt 0 ]; then
            paste -d '\\n' <(seq -f '>k%.0f' 1 "\$n_shared") work/shared.sorted > work/shared.fa
            printf '%s\\t%s\\n' "\$(stat -c %s work/shared.fa)" "\$PWD/work/shared.fa" > work/shared.list
            build_index_of_each_file work/shared.list 1 ${kmer_len} .w1.db
            shared_kmers=\$(db_kmers work/shared.fa.w1.db)
            subtract_db work/host.db work/shared.fa.w1.db ${db_out}
        else
            mv work/host.db ${db_out}
        fi
    else
        mv work/host.db ${db_out}
    fi
    final_kmers=\$(db_kmers ${db_out})
    rm -rf work

    printf '{\\n  "kmer_len": %d,\\n  "window": %d,\\n  "host_files": %d,\\n  "protect_files": %d,\\n  "host_kmers": %d,\\n  "shared_kmers_removed": %d,\\n  "kmers": %d\\n}\\n' \\
        ${kmer_len} ${window} ${[fasta].flatten().size()} ${protect_set.size()} \\
        "\$host_kmers" "\$shared_kmers" "\$final_kmers" > ${json_out}
    """

    stub:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def protect_set = protect ? [protect].flatten() : []
    """
    touch ${shellQuote("${prefix}.db")}
    printf '{\\n  "kmer_len": 32,\\n  "window": %d,\\n  "host_files": %d,\\n  "protect_files": %d,\\n  "host_kmers": 0,\\n  "shared_kmers_removed": 0,\\n  "kmers": 0\\n}\\n' \\
        ${window} ${[fasta].flatten().size()} ${protect_set.size()} > ${shellQuote("${prefix}.builddb.json")}
    """
}

//
// Single-quote a path for the shell, closing and reopening around any embedded
// quote, so a staged name or prefix containing whitespace stays one argument.
//
def shellQuote(path) {
    return "'" + path.toString().replace("'", "'\\''") + "'"
}
