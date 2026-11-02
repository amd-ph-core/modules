//
// The clockwork decontamination workflow for TB, simplifed
//   - minimap2 align to a database of known contaminants
//   - clockwork remove_contam removes read pairs that map to a contaminant
//     and regenerates fastq.gz files
//

include { CLOCKWORK_MINIMAP2     } from '../../../modules/amd-ph-core/clockwork/minimap2/main'
include { CLOCKWORK_REMOVECONTAM } from '../../../modules/amd-ph-core/clockwork/removecontam/main'

workflow FASTQ_CLOCKWORK_DECONTAMINATE {
    take:
    ch_reads         // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_reference     // path (mandatory): Path to reference genome for decontamination
    ch_ref_metadata  // path (mandatory): Path to reference metadata for decontamination

    main:
    // Intentionally empty, and kept only for interface compatibility. Both
    // modules in this subworkflow report through the `versions` TOPIC, which the
    // consuming pipeline collects globally via `channel.topic('versions')` --
    // topics are not scoped to a subworkflow, so nothing is lost by this channel
    // carrying nothing. Mixing them in here is not an option: the `emit: versions`
    // output they used to expose no longer exists.
    //
    // The emit stays because removing it would break any pipeline that reads
    // `FASTQ_CLOCKWORK_DECONTAMINATE.out.versions`. It can go once that is
    // confirmed unused.
    ch_versions = Channel.empty()

    CLOCKWORK_MINIMAP2 (
        ch_reads,
        ch_reference
    )

    CLOCKWORK_REMOVECONTAM (
        CLOCKWORK_MINIMAP2.out.sam,
        ch_ref_metadata
    )

    emit:
    counts               = CLOCKWORK_REMOVECONTAM.out.counts               // channel: [ val(meta), path(counts) ]
    decontaminated_fastq = CLOCKWORK_REMOVECONTAM.out.decontaminated_fastq // channel: [ val(meta), path(decontaminated_fastq) ]
    contaminant_fastq    = CLOCKWORK_REMOVECONTAM.out.contaminant_fastq    // channel: [ val(meta), path(contaminant_fastq) ]
    versions             = ch_versions                                     // channel: empty; versions arrive via the `versions` topic
}
