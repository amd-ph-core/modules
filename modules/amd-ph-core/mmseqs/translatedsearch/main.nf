// Tier 2 -- 6-frame translated search of the residual against a protein target DB.
//
// This is the highest-yield enrichment step measured: run on the residual of a 200,000-pair
// subsample it recruited 67.5% of what was left, at ~2,000 pairs/s and 2.4 GB peak on 14
// threads (99.8 s for 200,000 pairs). On the leftover pool it recovered 99.4-99.7% of what
// blastn -word_size 11 found and added a further 23-32%, with ZERO rRNA-only best hits across
// 3,000 adjudicated pairs even though rRNA ORFs were 90% of the target database by record
// count. That control is why the false-positive rate here is a measurement and not a hope.
//
// It exists as a separate tier because nucleotide alignment cannot reach these reads at all:
// at ~20% nucleotide divergence a 21-mer survives with probability 0.8^21 ~ 0.9%, so minimizer
// seeding never fires. Lowering the nucleotide identity floor from 80% to 50% recovered only
// 20-38% of them at flat ~77% precision, so this is not a threshold problem.
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
process MMSEQS_TRANSLATEDSEARCH {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // The tag this module carried in the binman pipeline --
    // `quay.io/biocontainers/mmseqs2:18.8cc5c--pl5321hd6d6fdc_0` -- does not exist. Docker
    // answers `manifest unknown`, so every real run of this module would have failed at pull
    // time; the measurements behind it were taken with a locally available mmseqs. The Seqera
    // community image is the same mmseqs 18.8cc5c, is what nf-core's own mmseqs modules use,
    // and sits on a registry the lint allow-list already permits.
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

    input:
    tuple val(meta), path(reads)
    path target_db

    output:
    tuple val(meta), path("${meta.id}.mmseqs.m8.gz"), emit: hits
    // See the note above the process about what may be written here.
    tuple val("${task.process}"), val('mmseqs2'), eval("mmseqs version"), emit: versions_mmseqs2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--search-type 2 --min-aln-len 30 -e 0.001'
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The target is passed as a prebuilt mmseqs DB rather than a FASTA: createdb on ~1.1M ORFs
    // is the dominant cost of a one-off search and repeating it per sample would swamp the
    // 100 s the search itself takes.
    def db_name = "${target_db}/targetdb"
    """
    set -euo pipefail

    # Mates are searched together: a hit on either mate places the fragment, and the assignment
    # step collapses mate names to fragment names anyway. Keeping them in one query DB avoids
    # paying the target index load twice.
    zcat ${reads} > query.fasta.tmp
    awk 'NR%4==1 {print ">" substr(\$0,2)} NR%4==2 {print}' query.fasta.tmp > query.fasta
    rm -f query.fasta.tmp

    mmseqs createdb query.fasta querydb

    mmseqs search \\
        querydb \\
        ${db_name} \\
        resdb \\
        tmp \\
        ${args} \\
        --threads ${task.cpus}

    mmseqs convertalis querydb ${db_name} resdb ${prefix}.mmseqs.m8 \\
        --format-output 'query,target,pident,alnlen,evalue,bits'

    gzip -n ${prefix}.mmseqs.m8
    rm -rf tmp resdb* querydb* query.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.mmseqs.m8.gz
    """
}
