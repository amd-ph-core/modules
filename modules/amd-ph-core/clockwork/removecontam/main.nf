process CLOCKWORK_REMOVECONTAM {
    tag "$meta.id"
    label 'process_low'

    container 'cdc-amd/crc:0.12.5_ed2e676_v10'

    input:
    tuple val(meta), path(contam_sam)
    path ref_metadata

    output:
    tuple val(meta), path("*.decontaminated.{1,2}.fastq.gz"), emit: decontaminated_fastq
    tuple val(meta), path("*.contaminants.{1,2}.fastq.gz")  , emit: contaminant_fastq
    tuple val(meta), path("*.counts.tsv")                   , emit: counts
    tuple val("${task.process}"), val('clockwork'), eval('clockwork version'), emit: versions_clockwork, topic: versions
    tuple val("${task.process}"), val('pysam'), eval('python -c "import pysam; print(pysam.__version__)"'), emit: versions_pysam, topic: versions
    tuple val("${task.process}"), val('pyfastaq'), eval('python -c "import pyfastaq; print(pyfastaq.__version__)"'), emit: versions_pyfastaq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def counts           = "${meta.id}.counts.tsv"
    def decontaminated_1 = "${meta.id}.decontaminated.1.fastq.gz"
    def decontaminated_2 = "${meta.id}.decontaminated.2.fastq.gz"
    def contaminants_1   = "${meta.id}.contaminants.1.fastq.gz"
    def contaminants_2   = "${meta.id}.contaminants.2.fastq.gz"

    """
    clockwork remove_contam \\
        $ref_metadata \\
        $contam_sam \\
        $counts \\
        $decontaminated_1 \\
        $decontaminated_2 \\
        --contam_out_1 $contaminants_1 \\
        --contam_out_2 $contaminants_2
    """

    stub:
    def counts           = "${meta.id}.counts.tsv"
    def decontaminated_1 = "${meta.id}.decontaminated.1.fastq.gz"
    def decontaminated_2 = "${meta.id}.decontaminated.2.fastq.gz"
    def contaminants_1   = "${meta.id}.contaminants.1.fastq.gz"
    def contaminants_2   = "${meta.id}.contaminants.2.fastq.gz"

    // The four FASTQs must be real gzip streams, not touched empty files. A
    // zero-byte `.gz` is not a valid member, so anything that opens it as one
    // throws -- nf-test's own md5 helper raises
    // `java.io.EOFException ... GZIPInputStream.readHeader` the moment a test
    // snapshots these outputs. Only `counts` is plain text and may be touched.
    //
    // `test_stub_gzip_syntax` did not catch this: it only inspects lines whose
    // last characters are a literal `.gz`, and the old `touch` line ended in a
    // Groovy variable, so the rule skipped it entirely.
    """
    touch $counts
    echo "" | gzip > $decontaminated_1
    echo "" | gzip > $decontaminated_2
    echo "" | gzip > $contaminants_1
    echo "" | gzip > $contaminants_2
    """
}
