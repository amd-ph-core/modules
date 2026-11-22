process SRAHUMANSCRUBBER_BUILDDB {
    tag "${meta.id}"
    label 'process_high'

    // The STAT database tools (build_index_of_each_file, merge_db, subtract_db,
    // aligns_to) come from ncbi/ngs-tools (tools/tax) and have no conda package or
    // BioContainer: the bioconda sra-human-scrubber package ships aligns_to only.
    // The AMDP sra-human-scrubber image carries them alongside the scrubber, built
    // from the same aligns_to, so one image serves this module and scrub.
    // Consumers override the registry and the exact tag via their own config, as
    // with the other AMDP images.
    container 'oamd-bio-sra-human-scrubber:2.2.1'

    input:
    tuple val(meta), path(fasta, stageAs: "host/*")
    path(protect, stageAs: "protect/*")
    val window

    output:
    tuple val(meta), path("*.db")                    , emit: db
    tuple val(meta), path("*.builddb.json")          , emit: json
    tuple val("${task.process}"), val('aligns_to'), eval('aligns_to 2>&1 | sed -n "s/.*aligns_to version //p" | head -n 1 | grep . || echo unknown'), topic: versions, emit: versions_aligns_to
    tuple val("${task.process}"), val('build_index_of_each_file'), eval('build_index_of_each_file 2>&1 | sed -n "s/.*build_index_of_each_file version //p" | head -n 1 | grep . || echo unknown'), topic: versions, emit: versions_build_index_of_each_file

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def protect_set = protect ? [protect].flatten() : []
    def kmer_len    = 32
    """
    # A STAT .db is a header (format version, k) and a sorted set of canonical
    # k-mers. From each consecutive window of a sequence, build_index_of_each_file
    # keeps the k-mer with the smallest hash: window 1 keeps every k-mer, window 32
    # one in 32 (NCBI's human_filter.db is consistent with 32). aligns_to tests
    # every k-mer of a read, so an error-free host read of at least 2 x window +
    # k - 1 bases always contains a kept k-mer.
    db_kmers() { echo \$(( (\$(stat -c %s "\$1") - 24) / 8 )); }

    # One plain FASTA per input; the tools do not read gzip.
    mkdir -p plain
    for f in host/*; do
        gzip -cdf "\$f" > "plain/\$(basename "\${f%.gz}")"
    done
    : > host.list
    for f in plain/*; do
        printf '%s\\t%s\\n' "\$(stat -c %s "\$f")" "\$PWD/\$f" >> host.list
    done
    build_index_of_each_file host.list ${window} ${kmer_len} .w${window}.db
    # merge_db sizes its hash table up front: about two slots per expected k-mer.
    expected=\$(awk -v w=${window} '{ s += \$1 } END { printf "%d", 2 * s / w + 1000 }' host.list)
    merge_db host.list .w${window}.db host.db "\$expected"
    host_kmers=\$(db_kmers host.db)

    shared_kmers=0
    if [ -n "${protect_set.join(' ')}" ]; then
        # Protecting sequences: every host k-mer found in them is dropped, so a read
        # from them cannot be identified as host through a shared k-mer. aligns_to
        # streams the protected FASTAs against the host set and prints each host
        # k-mer it meets; that list is rebuilt into a set at window 1 (which also
        # puts each k-mer in canonical orientation) and subtracted.
        : > shared.txt
        for f in ${protect_set.join(' ')}; do
            gzip -cdf "\$f" > protect.fa
            aligns_to -db host.db -num_threads ${task.cpus} -print_kmers_only protect.fa >> shared.txt
        done
        rm -f protect.fa
        sort -u shared.txt | awk '{ print ">k" NR; print }' > shared.fa
        rm shared.txt
        if [ -s shared.fa ]; then
            printf '%s\\t%s\\n' "\$(stat -c %s shared.fa)" "\$PWD/shared.fa" > shared.list
            build_index_of_each_file shared.list 1 ${kmer_len} .w1.db
            shared_kmers=\$(db_kmers shared.fa.w1.db)
            subtract_db host.db shared.fa.w1.db ${prefix}.db
        else
            mv host.db ${prefix}.db
        fi
    else
        mv host.db ${prefix}.db
    fi
    final_kmers=\$(db_kmers ${prefix}.db)
    rm -rf plain host.list host/*.w${window}.db plain/*.w${window}.db shared.* host.db

    printf '{\\n  "kmer_len": %d,\\n  "window": %d,\\n  "host_files": %d,\\n  "protect_files": %d,\\n  "host_kmers": %d,\\n  "shared_kmers_removed": %d,\\n  "kmers": %d\\n}\\n' \\
        ${kmer_len} ${window} ${[fasta].flatten().size()} ${protect_set.size()} \\
        "\$host_kmers" "\$shared_kmers" "\$final_kmers" > ${prefix}.builddb.json
    """

    stub:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def protect_set = protect ? [protect].flatten() : []
    """
    touch ${prefix}.db
    printf '{\\n  "kmer_len": 32,\\n  "window": %d,\\n  "host_files": %d,\\n  "protect_files": %d,\\n  "host_kmers": 0,\\n  "shared_kmers_removed": 0,\\n  "kmers": 0\\n}\\n' \\
        ${window} ${[fasta].flatten().size()} ${protect_set.size()} > ${prefix}.builddb.json
    """
}
