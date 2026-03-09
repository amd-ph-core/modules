include { BWA_INDEX                                      } from '../../../modules/nf-core/bwa/index/main'
include { PICARD_CREATESEQUENCEDICTIONARY                } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                                 } from '../../../modules/nf-core/samtools/faidx/main'
include { VARPIPE_BWAMEM                                 } from '../../../modules/amd-ph-core/varpipe/bwamem/main'
include { PICARD_SAMFORMATCONVERTER                      } from '../../../modules/amd-ph-core/picard/samformatconverter/main'
include { PICARD_SORTSAM                                 } from '../../../modules/amd-ph-core/picard/sortsam/main'
include { PICARD_MARKDUPLICATES                          } from '../../../modules/nf-core/picard/markduplicates/main'
include { PICARD_BUILDBAMINDEX                           } from '../../../modules/amd-ph-core/picard/buildbamindex/main'
include { PICARD_BUILDBAMINDEX as PICARD_BUILDBAMINDEX_2 } from '../../../modules/amd-ph-core/picard/buildbamindex/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_UNMAPPED        } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW                                  } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_MAPPED          } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_DEPTH                                 } from '../../../modules/nf-core/samtools/depth/main'

workflow VARPIPE_PROCESSBAM {
    take:
    ch_reads     // channel (mandatory): [ val(meta), path(reads) ]
    ch_ref_fasta // channel (mandatory): [ val(meta), path(ref_fasta) ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_inputs = Channel.empty()

    // ALIGNERS
    BWA_INDEX (
        ch_ref_fasta
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    PICARD_CREATESEQUENCEDICTIONARY (
        ch_ref_fasta
    )
    ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)

    SAMTOOLS_FAIDX (
        ch_ref_fasta,
        [[],[]],
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    ch_index = BWA_INDEX.out.index
        .join(PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict)
        .join(SAMTOOLS_FAIDX.out.fai)
        .join(ch_ref_fasta)
        .collect()

    VARPIPE_BWAMEM (
        ch_reads,
        ch_index
    )
    ch_versions = ch_versions.mix(VARPIPE_BWAMEM.out.versions)

    // CONVERT SAM TO BAM
    PICARD_SAMFORMATCONVERTER (
        VARPIPE_BWAMEM.out.sam
    )
    ch_versions = ch_versions.mix(PICARD_SAMFORMATCONVERTER.out.versions)

    // RUN MAPPING REPORT AND MARK DUPLICATES
    PICARD_SORTSAM (
        PICARD_SAMFORMATCONVERTER.out.bam,
        'coordinate'
    )
    ch_versions = ch_versions.mix(PICARD_SORTSAM.out.versions)

    PICARD_MARKDUPLICATES (
        PICARD_SORTSAM.out.bam,
        ch_ref_fasta,
        SAMTOOLS_FAIDX.out.fai.collect()
    )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
    ch_multiqc_inputs = ch_multiqc_inputs.mix(PICARD_MARKDUPLICATES.out.metrics
        .map { meta, metrics -> metrics })

    PICARD_BUILDBAMINDEX (
        PICARD_MARKDUPLICATES.out.bam
    )
    ch_versions = ch_versions.mix(PICARD_BUILDBAMINDEX.out.versions)

    SAMTOOLS_VIEW_UNMAPPED (
        PICARD_MARKDUPLICATES.out.bam.join(PICARD_BUILDBAMINDEX.out.bai),
        [[],[]],
        [],
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_UNMAPPED.out.versions)

    SAMTOOLS_VIEW (
        PICARD_MARKDUPLICATES.out.bam.join(PICARD_BUILDBAMINDEX.out.bai),
        [[],[]],
        [],
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    PICARD_BUILDBAMINDEX_2 (
        SAMTOOLS_VIEW.out.bam
    )
    ch_versions = ch_versions.mix(PICARD_BUILDBAMINDEX_2.out.versions)

    SAMTOOLS_VIEW_MAPPED (
        SAMTOOLS_VIEW.out.bam.join(PICARD_BUILDBAMINDEX_2.out.bai),
        [[],[]],
        [],
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_MAPPED.out.versions)

    // RUN GENOME COVERAGE STATISTICS
    SAMTOOLS_DEPTH (
        SAMTOOLS_VIEW.out.bam,
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions)

    emit:
    bam                   = SAMTOOLS_VIEW.out.bam                              // channel: [ val(meta), path(bam) ]
    bai                   = PICARD_BUILDBAMINDEX_2.out.bai                     // channel: [ val(meta), path(bai) ]
    coverage              = SAMTOOLS_DEPTH.out.tsv                             // channel: [ val(meta), path(tsv) ]
    fai                   = SAMTOOLS_FAIDX.out.fai                             // channel: [ val(meta), path(fai) ]
    unmapped              = SAMTOOLS_VIEW_UNMAPPED.out.bam                     // channel: [ val(meta), path(bam) ]
    mapped                = SAMTOOLS_VIEW_MAPPED.out.bam                       // channel: [ val(meta), path(bam) ]
    reference_dict        = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict // channel: [ val(meta), path(ref_dict) ]
    multiqc_inputs        = ch_multiqc_inputs                                  // channel: [ path(markduplicates_metrics) ]
    versions              = ch_versions                                        // channel: [ versions.yml ]
}
