process PICARD_SAMFORMATCONVERTER {

    tag "${meta.id}"
    label 'process_medium'

    container 'cdc-amd/picard:3.4.0_ac49616_v8'

    input:
    tuple val(meta), path(alnSam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('picard'), eval('picard SamFormatConverter --version 2>&1 | sed -n "s/^Version://p"'), emit: versions_picard, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The output name comes from `ext.prefix` alone. It used to take a second
    // key, `ext.output_suffix`, which nf-core lint rejects outright:
    // `main_nf_ext_key` permits only args/args<N>/prefix/prefix2/use_gpu, and a
    // module-repo lint rule cannot be waived, so the key had to go rather than
    // be silenced.
    //
    // Folding it into `ext.prefix` also fixes the name it produced. The one
    // caller that set it passed `_GATK.bam` into `OUTPUT=${prefix}.${suffix}`,
    // so the file came out as `<id>._GATK.bam` -- a stray dot ahead of the
    // underscore. Carrying the discriminator in `ext.prefix` gives
    // `<id>_GATK.bam`, which is the intended name and matches what every other
    // process in that same config already does (`ext.prefix = {
    // "${meta.id}_GATK_sdr" }`). Downstream consumers glob `*.bam`, so the
    // rename does not change what they match.
    """
    picard SamFormatConverter \\
        INPUT=${alnSam} \\
        VALIDATION_STRINGENCY=${args} \\
        OUTPUT=${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
