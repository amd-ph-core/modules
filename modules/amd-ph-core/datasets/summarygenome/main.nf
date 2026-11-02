// Enumerate every assembly under a parent taxon in a single `datasets summary`
// call. No genome data is fetched here, so this stays cheap enough to run up
// front and let downstream filtering decide what is worth downloading.
//
// Ported from securebio/nao-mgs-workflow ENUMERATE_VIRAL_ACCESSIONS, minus its
// `dataformat tsv genome` flattening step. `dataformat` has no public source --
// it ships only as an NCBI prebuilt -- so an image that builds `datasets` from
// source does not carry it, which is the case for the OAMD ECR build even though
// both the StaPH-B image and the conda package do bundle it. Staying off
// `dataformat` is what keeps this module runnable on all three. The raw JSON
// Lines report carries the same fields, so flattening belongs downstream in a
// process that has a JSON parser.
process DATASETS_SUMMARYGENOME {

    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Fully qualified on purpose. The bare `staphb/ncbi-datasets` form resolved
    // against whatever `docker.registry` happened to be configured -- the repo
    // test config sets it to a mirror -- so the same tag could pull a different
    // image depending on the consumer's config, and it forced `staphb` onto the
    // lint `--registry` allow-list as an unqualified prefix. StaPH-B publishes
    // this image to quay.io, which is nf-core's canonical registry and already
    // permitted, so naming the registry costs nothing and pins what actually runs.
    container 'quay.io/staphb/ncbi-datasets:18.35.0'

    input:
    tuple val(meta), val(taxon)

    output:
    tuple val(meta), path("*.assembly_report.jsonl.gz"), emit: report
    tuple val(meta), path("*.accessions.txt")          , emit: accessions
    tuple val("${task.process}"), val('ncbi-datasets-cli'), eval('datasets --version 2>&1 | sed "s/datasets version: //"'), emit: versions_ncbidatasetscli, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--assembly-source all'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail

    # 1. Enumerate all assemblies under the taxon.
    datasets summary genome taxon ${taxon} \\
        --as-json-lines \\
        ${args} \\
        > ${prefix}.assembly_report.jsonl

    # 2. Pull the accession list out of the report. `accession` is the first key
    # of every record, so first-match per line takes the assembly's own
    # accession and not the nested `paired_accession` / `current_accession`.
    # 13 = length of '"accession":"'; the extra char dropped is its closing quote.
    awk 'match(\$0, /"accession":"[^"]*"/) { print substr(\$0, RSTART+13, RLENGTH-14) }' \\
        ${prefix}.assembly_report.jsonl \\
        | LC_ALL=C sort \\
        > ${prefix}.accessions.txt

    # A taxon with no assemblies is a bad input, not an empty result set. Fail
    # here rather than let downstream build an empty genome DB.
    if [ ! -s ${prefix}.accessions.txt ]; then
        echo "No assemblies found for taxon ${taxon}" >&2
        exit 1
    fi

    gzip -n ${prefix}.assembly_report.jsonl
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // `echo "" | gzip >` is the only form the linter accepts -- its regex is
    // `echo\s+""\s*\|\s*gzip\s*>\s*.*\.gz$`, which rejects the `-n` flag even
    // though `-n` is what the real script uses to keep the output reproducible.
    """
    echo "" | gzip > ${prefix}.assembly_report.jsonl.gz
    touch ${prefix}.accessions.txt
    """
}
