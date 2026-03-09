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
    path "versions.yml"                                     , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clockwork: \$(clockwork version)
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
        pyfastaq: \$(python -c "import pyfastaq; print(pyfastaq.__version__)")
    END_VERSIONS
    """

    stub:
    def counts           = "${meta.id}.counts.tsv"
    def decontaminated_1 = "${meta.id}.decontaminated.1.fastq.gz"
    def decontaminated_2 = "${meta.id}.decontaminated.2.fastq.gz"
    def contaminants_1   = "${meta.id}.contaminants.1.fastq.gz"
    def contaminants_2   = "${meta.id}.contaminants.2.fastq.gz"

    """
    touch $counts $decontaminated_1 $decontaminated_2 $contaminants_1 $contaminants_2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clockwork: \$(clockwork version)
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
        pyfastaq: \$(python -c "import pyfastaq; print(pyfastaq.__version__)")
    END_VERSIONS
    """
}
