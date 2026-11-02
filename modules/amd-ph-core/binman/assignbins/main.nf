// Fold tier-2 and tier-3 enrichment hits into the tier-1 kraken2 family bins.
//
// The policy and its measured justification live in the template script. In short: union the
// family labels a read earned from any tier, place it if exactly one family survives, hold it in
// the residual if more than one does. Enrichment never overrides a tier-1 bin, because no tier
// is a superset of any other -- kraken2 alone held 1,706 reads per 200,000 that nothing else
// recovered.
process BINMAN_ASSIGNBINS {

    tag "${meta.id}"
    label 'process_medium'

    // Public, and deliberately not the AMDP `oamd-bio-python:3.12-orb` image this module used
    // while it lived in the binman pipeline. The script is pure standard library -- argparse,
    // collections, csv, gzip, pathlib, shlex -- so it has no claim on a bespoke image, and a
    // module offered to other pipelines should not drag an ECR dependency along with it. The
    // tag is a minor version rather than a build hash, because that is the most specific thing
    // biocontainers publishes for this image; the patch version it reports is free to move
    // within 3.12, which is why the test snapshot does not pin it.
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(residual), path(bin_manifest), path(mmseqs_hits), path(rammap_paf)
    path rammap_blocks

    output:
    tuple val(meta), path("enriched/*.fastq.gz")             , emit: bins    , optional: true
    tuple val(meta), path("${meta.id}.residual_*.fastq.gz")  , emit: residual
    tuple val(meta), path("${meta.id}.assignment_report.tsv"), emit: report
    // Declared as the `versions` TOPIC while still being a real file. This module's script is a
    // Python `template`, and Nextflow rejects an `eval` output on any non-Bash interpreter --
    // "Process output of type 'eval' is only allowed with Bash process scripts". Lint accepts
    // `path "versions.yml"` on a topic line for exactly this case, so the template writes the
    // file and the topic carries it.
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Thresholds and the cross-mapping policy reach the script through `task.ext.args`, which
    // the template splits and appends to argv. They used to be read from `params` directly;
    // that cannot survive promotion to a shared repo, because a consumer pipeline that has
    // never declared `params.rammap_min_identity` would fail when the script is evaluated.
    // The script's own argparse defaults -- 60 bp, 80% identity, hold cross-mapped reads --
    // are the measured settings, so a consumer that sets no args gets the validated behaviour.
    template 'assign_bins.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mate_list = meta.single_end ? '1' : '1 2'
    """
    mkdir -p enriched
    echo "" | gzip > enriched/taxid_12058_picornaviridae_enriched_R1.fastq.gz
    for m in ${mate_list}; do
        echo "" | gzip > ${prefix}.residual_\$m.fastq.gz
    done
    printf '%s\\n' 'metric\tvalue' 'enriched_into_bins\t0' > ${prefix}.assignment_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
