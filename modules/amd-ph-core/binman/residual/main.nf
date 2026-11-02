// Build the residual pool that tiers 2 and 3 search.
//
// The residual is two streams because kraken2-hunter emits it as two, and that split is a
// property of the patch rather than a choice: its bin writer sits inside `if (call)`, so a read
// kraken2 could not classify never reaches the writer and leaves via --unclassified-out, while
// a read whose LCA sits above the bin rank does reach it and lands in the unresolved bin. Both
// are "classified into no family bin", which is what the enrichment tiers operate on.
//
// Measured proportions on an rRNA-depleted wastewater library: the unclassified stream is 72.4%
// of the depleted library and the unresolved bin is 0.19%. The LCA bin being nearly empty is
// itself the finding that made "kraken2 family bin plus LCA bin" indistinguishable from "kraken2
// family bin" -- but those reads are still placeable, so they are kept.
//
// This module does NOT glob for its inputs, and that is a correction rather than a style
// preference. An earlier version assembled each mate with a shell glob, and the two producers
// do not share a naming convention: the bin writer emits `unresolved_family_R1.fastq.gz` for
// paired input but `unresolved_family.fastq.gz` for single, while kraken2 substitutes `#` with
// `_1`/`_2` rather than `1`/`2`. Two of the four combinations missed, so the pool was built from
// one stream instead of two and still looked plausible -- 182 reads for a sample whose
// unresolved bin held exactly 182. Using the staged inputs directly removes the guesswork, and
// the read-count assertion below turns any remaining mismatch into a failure instead of a
// quietly truncated pool.
// A note on the comment inside the declared-outputs section below, because it has now cost two
// rounds of debugging: nf-core's linter parses COMMENT lines in that section as if they were
// declarations, so a comment there must not contain any of the block keywords. Writing the
// channel-topic keyword in backticks registered a phantom topic named `versions`` -- trailing
// backtick included -- which failed correct_meta_topics, and `--fix` then wrote that broken name
// into meta.yml as a real key rather than reporting it. Replacing it with prose that happened to
// mention the outputs keyword was worse: the parser took it as the start of a new section and
// found zero topics and one output. Explanations therefore live out here, where nothing parses
// them, and the comment in there says nothing.
//
// The version is reported on the versions channel topic, which is what nf-core 4.x collects. A
// legacy versions.yml emit alongside it would trip main_nf_version_topic, which compares counts.
process BINMAN_RESIDUAL {

    tag "${meta.id}"
    label 'process_low'

    // Public, and deliberately not the `cdc-amd/ubuntu` mirror this module used while it lived
    // in the binman pipeline. All it needs is cat, gzip and wc, so a module offered to other
    // pipelines has no reason to require access to our registry.
    container 'quay.io/nf-core/ubuntu:24.04'

    input:
    tuple val(meta), path(unresolved), path(unclassified)

    output:
    tuple val(meta), path("${meta.id}.residual_pool_*.fastq.gz"), emit: residual
    // See the note above the process about what may be written here.
    tuple val("${task.process}"), val('gzip'), eval("gzip --version | sed -n '1s/.* //p'"), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def n_mates = meta.single_end ? 1 : 2
    // Each group arrives as [] when absent, a bare Path when it staged one file, or a List when
    // it staged several. Normalising through [x].flatten() is required, not defensive: a bare
    // java.nio.file.Path is Iterable over its own name SEGMENTS, so calling .sort{} on one
    // returns the path components rather than a one-element list.
    def unres = [unresolved].flatten().findAll { it }.sort { it.name }
    def unc = [unclassified].flatten().findAll { it }.sort { it.name }
    // A group may legitimately be absent -- a library with an empty unresolved bin is normal.
    // A group that is present but has the wrong number of files is not, and must stop the run.
    if (unres && unres.size() != n_mates) {
        error "BINMAN_RESIDUAL(${meta.id}): unresolved bin has ${unres.size()} file(s), " +
              "expected ${n_mates} for single_end=${meta.single_end}: ${unres*.name}"
    }
    if (unc && unc.size() != n_mates) {
        error "BINMAN_RESIDUAL(${meta.id}): unclassified stream has ${unc.size()} file(s), " +
              "expected ${n_mates} for single_end=${meta.single_end}: ${unc*.name}"
    }
    if (!unres && !unc) {
        error "BINMAN_RESIDUAL(${meta.id}): neither an unresolved bin nor an unclassified " +
              "stream was produced. On a real library the unclassified stream alone is " +
              "typically 60-80% of reads, so this means the upstream call failed, not that " +
              "the residual is genuinely empty."
    }
    // Sorted by name, each group's files are in mate order, so index i is mate i+1.
    def parts = (0..<n_mates).collect { i ->
        def srcs = []
        if (unres) { srcs << unres[i] }
        if (unc) { srcs << unc[i] }
        "cat ${srcs.join(' ')} > ${prefix}.residual_pool_${i + 1}.fastq.gz"
    }.join('\n    ')
    def inputs_m1 = ((unres ? [unres[0]] : []) + (unc ? [unc[0]] : [])).join(' ')
    """
    set -euo pipefail

    ${parts}

    # Concatenated gzip members form a valid gzip stream, so no recompression is needed. Verify
    # that rather than assume it.
    for f in ${prefix}.residual_pool_*.fastq.gz; do
        gzip -t "\$f"
    done

    # Account for every input read. A glob or a name mismatch that silently contributes nothing
    # is the failure mode this module already had once, and it produced numbers that looked
    # entirely reasonable. Comparing mate-1 line counts makes that unmissable.
    want=\$(zcat ${inputs_m1} | wc -l)
    got=\$(zcat ${prefix}.residual_pool_1.fastq.gz | wc -l)
    if [ "\$want" -ne "\$got" ]; then
        echo "BINMAN_RESIDUAL(${meta.id}): mate 1 residual pool has \$got lines but its inputs have \$want" >&2
        exit 1
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mate_list = meta.single_end ? '1' : '1 2'
    """
    for m in ${mate_list}; do
        echo "" | gzip > ${prefix}.residual_pool_\$m.fastq.gz
    done
    """
}
