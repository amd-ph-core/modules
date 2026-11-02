process FADSO {
    tag "${meta.id}"
    label 'process_low'

    // Upstream (hkuwahara/sleaping) publishes no releases or tags, so the
    // container is built from a pinned commit and there is no conda package.
    // No public image exists: fadso has no conda package and no BioContainer,
    // and upstream publishes no releases. This is the AMDP ECR repository being
    // created for it; consumers override the registry via their own config.
    container 'oamd-bio-fadso:20230810_1fedb59_v0'

    input:
    tuple val(meta), path(reads)
    val max_reads
    val seed

    output:
    tuple val(meta), path("*.downsampled.fastq.gz"), emit: reads
    // fadso has no --version flag and upstream publishes no releases, so the
    // pinned source commit stands in for a version per ph-core guidance. It is
    // echoed rather than hardcoded as a val so the output matches the shape
    // nf-core expects for the versions topic.
    tuple val("${task.process}"), val('fadso'), eval("echo 1fedb59"), topic: versions, emit: versions_fadso

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        fadso single \\
            ${args} \\
            -r ${seed} \\
            -k ${max_reads} \\
            -i ${reads} \\
            -o ${prefix}.downsampled.fastq.gz \\
            -z
        """
    }
    else {
        """
        fq=( ${reads} )
        fadso pair \\
            ${args} \\
            -r ${seed} \\
            -k ${max_reads} \\
            -1 "\${fq[0]}" \\
            -2 "\${fq[1]}" \\
            -a ${prefix}_1.downsampled.fastq.gz \\
            -b ${prefix}_2.downsampled.fastq.gz \\
            -z
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outputs = meta.single_end
        ? "${prefix}.downsampled.fastq.gz"
        : "${prefix}_1.downsampled.fastq.gz ${prefix}_2.downsampled.fastq.gz"
    // Each output is gzipped directly rather than copied from a shared temp file.
    // The previous form was functionally fine but tripped `test_stub_gzip_syntax`:
    // that rule inspects every line ending in `.gz` and demands it match
    // `echo "" | gzip > ...`, and `rm tmp.gz` ends in `.gz` without matching. Writing
    // each file in the loop removes the temp file and the false positive together.
    """
    for out in ${outputs}; do
        echo "" | gzip > "\$out"
    done
    """
}
