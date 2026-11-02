// Classify reads and write per-rank taxonomic bins in the SAME pass.
//
// This is Hunter Seabolt's fork of kraken2 (hseabolt/kraken2 @ 23af0bd, on upstream 2.17.1),
// which adds an output-layer extension: as classify walks the reads it appends each fragment
// to a bin chosen by the nearest ancestor at a requested rank. Nothing about the database
// format or the classification call changes -- the same reads get the same taxids -- so the
// report this emits is comparable with plain kraken2 output.
//
// Why it is worth a separate module rather than kraken2 + a splitting script: the alternative
// is re-reading the FASTQs and the per-read output to partition them afterwards, which on a
// metagenome means a second full pass over the reads for no new information. Measured here at
// 100,000 pairs against the 16 GB standard database, binning to family cost 7.6 s wall and
// produced 16 family bins plus an unresolved bin.
//
// CONDA IS DELIBERATELY ABSENT. There is no conda package for this fork, and declaring
// `conda "${moduleDir}/environment.yml"` would resolve to upstream bioconda::kraken2 -- which
// classifies correctly but has no --bin-out flag at all. The process would then fail with an
// unrecognised-option error that looks like a module bug rather than a missing patch, or worse,
// a wrapper could ignore it and silently produce no bins. A container-only module fails
// honestly instead.
process KRAKEN2HUNTER_CLASSIFY {

    tag "${meta.id}"
    label 'process_high'

    container '025066257930.dkr.ecr.us-east-1.amazonaws.com/oamd-bio-kraken2-hunter:2.17.1-23af0bd'

    input:
    tuple val(meta), path(reads)
    path db
    val bin_rank

    output:
    tuple val(meta), path("${meta.id}.kraken2.report.txt")          , emit: report
    tuple val(meta), path("bins/bins.tsv")                          , emit: bin_manifest
    tuple val(meta), path("bins/taxid_*.fastq.gz")                  , emit: bins              , optional: true
    tuple val(meta), path("bins/unresolved_*.fastq.gz")             , emit: unresolved        , optional: true
    tuple val(meta), path("${meta.id}.unclassified*.fastq.gz")      , emit: unclassified      , optional: true
    tuple val(meta), path("${meta.id}.kraken2.output.txt.gz")       , emit: classified_reads  , optional: true
    // Version is captured by `eval` at task time rather than written to a file. `sed -n s//p`
    // strips the "Kraken version " banner prefix, leaving a bare 2.17.1 with no leading "v" --
    // ph-core requires the version number alone.
    tuple val("${task.process}"), val('kraken2-hunter'), eval('kraken2 --version | sed -n "s/^Kraken version //p"'), emit: versions_kraken2hunter, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? '' : '--paired'
    // kraken2 substitutes '#' with '_1'/'_2' -- WITH the underscore, not a bare mate number.
    // Writing `unclassified_#` therefore yields `unclassified__1`, a double underscore that a
    // downstream name match will miss. Verified against the built container rather than assumed.
    def unc = meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    """
    set -euo pipefail

    mkdir -p bins

    kraken2 \\
        --db ${db} \\
        --threads ${task.cpus} \\
        ${paired} \\
        --gzip-compressed \\
        --report ${prefix}.kraken2.report.txt \\
        --output ${prefix}.kraken2.output.txt \\
        --unclassified-out '${unc}' \\
        --bin-out bins \\
        --bin-rank ${bin_rank} \\
        ${args} \\
        ${reads}

    # The bins themselves are written compressed by the tool (--bin-compression gzip, supplied
    # through ext.args). These two streams are not, so compress them here rather than leave an
    # uncompressed metagenome-sized FASTQ in the work directory.
    gzip -n ${prefix}.kraken2.output.txt
    for f in ${prefix}.unclassified*.fastq; do
        [ -e "\$f" ] || continue
        gzip -n "\$f"
    done

    # bins.tsv is the tool's own manifest: bin taxid, rank, scientific name, fragment/read/base
    # counts and the file for each mate. Downstream joins on it instead of parsing filenames,
    # because the taxid is authoritative and the sanitised scientific name in the filename is
    # not round-trippable.
    test -s bins/bins.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p bins
    touch ${prefix}.kraken2.report.txt
    printf '%s\\n' 'bin_taxid\trank\tscientific_name\tfragments\treads\tbases\tsequence_file_1\tsequence_file_2' > bins/bins.tsv
    echo "" | gzip > bins/taxid_12058_Picornaviridae_R1.fastq.gz
    echo "" | gzip > bins/taxid_12058_Picornaviridae_R2.fastq.gz
    echo "" | gzip > bins/unresolved_family_R1.fastq.gz
    echo "" | gzip > bins/unresolved_family_R2.fastq.gz
    echo "" | gzip > ${prefix}.unclassified_1.fastq.gz
    echo "" | gzip > ${prefix}.unclassified_2.fastq.gz
    echo "" | gzip > ${prefix}.kraken2.output.txt.gz
    """
}
