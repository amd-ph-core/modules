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
    ch_versions = Channel.empty()

    CLOCKWORK_MINIMAP2 (
        ch_reads,
        ch_reference
    )
    ch_versions = ch_versions.mix(CLOCKWORK_MINIMAP2.out.versions.first())

    CLOCKWORK_REMOVECONTAM (
        CLOCKWORK_MINIMAP2.out.sam,
        ch_ref_metadata
    )
    ch_versions = ch_versions.mix(CLOCKWORK_REMOVECONTAM.out.versions.first())

    emit:
    counts               = CLOCKWORK_REMOVECONTAM.out.counts               // channel: [ val(meta), path(counts) ]
    decontaminated_fastq = CLOCKWORK_REMOVECONTAM.out.decontaminated_fastq // channel: [ val(meta), path(decontaminated_fastq) ]
    contaminant_fastq    = CLOCKWORK_REMOVECONTAM.out.contaminant_fastq    // channel: [ val(meta), path(contaminant_fastq) ]
    versions             = ch_versions                                     // channel: [ versions.yml ]
}
