include { GATK4_MUTECT2 as GATK4_MUTECT2_ALL                                        } from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_LEFTALIGNANDTRIMVARIANTS as GATK4_LEFTALIGNANDTRIMVARIANTS_ALL      } from '../../../modules/nf-core/gatk4/leftalignandtrimvariants/main'
include { GATK4_FILTERMUTECTCALLS as GATK4_FILTERMUTECTCALLS_ALL                    } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'
include { SNPEFF_SNPEFF as VARPIPE_SNPEFF_ALL                                       } from '../../../modules/nf-core/snpeff/snpeff/main'

include { GATK4_MUTECT2 as GATK4_MUTECT2_TARGETED                                   } from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_LEFTALIGNANDTRIMVARIANTS as GATK4_LEFTALIGNANDTRIMVARIANTS_TARGETED } from '../../../modules/nf-core/gatk4/leftalignandtrimvariants/main'
include { GATK4_FILTERMUTECTCALLS as GATK4_FILTERMUTECTCALLS_TARGETED               } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'
include { SNPEFF_SNPEFF as VARPIPE_SNPEFF_TARGETED                                  } from '../../../modules/nf-core/snpeff/snpeff/main'

workflow VARPIPE_VARIANTANALYSIS {
    take:
    ch_input_all        // channel (mandatory): [ val(meta), path(bam), path(bai), [] ]
    ch_input_targeted   // channel (mandatory): [ val(meta), path(bam), path(bai), [] ]
    ch_fai              // channel (mandatory): [ val(meta), path(fai) ]
    ch_ref_dict         // channel (mandatory): [ val(meta), path(reference_dict) ]
    ch_ref_fasta        // channel (mandatory): [ val(meta),path(ref_fasta) ]
    snpeff_db           // string (mandatory): SnpEff database name (string) for variant annotation
    snpeff_cache        // path (mandatory): Path to SnpEff cache directory containing annotation databases

    main:
    ch_versions = Channel.empty()
    ch_multiqc_inputs = Channel.empty()

    // FULL GENOME ANALYSIS
    GATK4_MUTECT2_ALL (
        ch_input_all,
        ch_ref_fasta,
        ch_fai,
        ch_ref_dict,
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GATK4_MUTECT2_ALL.out.versions)

    ch_leftalignandtrimvariants_all = GATK4_MUTECT2_ALL.out.vcf
        .join(GATK4_MUTECT2_ALL.out.tbi)
        .map { meta, vcf, tbi ->
            return [ meta, vcf, tbi, [] ]
        }

    GATK4_LEFTALIGNANDTRIMVARIANTS_ALL (
        ch_leftalignandtrimvariants_all,
        ch_ref_fasta.map{ meta, fasta -> [ fasta ] },
        ch_fai.map{ meta, fai -> [ fai ] },
        ch_ref_dict.map{ meta, dict -> [ dict ] }
    )
    ch_versions = ch_versions.mix(GATK4_LEFTALIGNANDTRIMVARIANTS_ALL.out.versions)

    ch_filtermutectcalls_all = GATK4_LEFTALIGNANDTRIMVARIANTS_ALL.out.vcf
        .join(GATK4_LEFTALIGNANDTRIMVARIANTS_ALL.out.tbi)
        .join(GATK4_MUTECT2_ALL.out.stats)
        .map { meta, vcf, tbi, stats ->
            return tuple(
                meta,
                vcf,
                tbi,
                stats,
                [],
                [],
                [],
                []
            )
        }

    GATK4_FILTERMUTECTCALLS_ALL (
        ch_filtermutectcalls_all,
        ch_ref_fasta,
        ch_fai,
        ch_ref_dict
    )
    ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS_ALL.out.versions)

    VARPIPE_SNPEFF_ALL (
        GATK4_FILTERMUTECTCALLS_ALL.out.vcf,
        snpeff_db,
        snpeff_cache.map { cache -> [[], cache] }
    )
    ch_versions = ch_versions.mix(VARPIPE_SNPEFF_ALL.out.versions)
    ch_multiqc_inputs = ch_multiqc_inputs.mix(VARPIPE_SNPEFF_ALL.out.report.map { meta, report -> report })

    // DR LOCI ANALYSIS
    GATK4_MUTECT2_TARGETED (
        ch_input_targeted,
        ch_ref_fasta,
        ch_fai,
        ch_ref_dict,
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GATK4_MUTECT2_TARGETED.out.versions)

    ch_leftalignandtrimvariants_targeted = GATK4_MUTECT2_TARGETED.out.vcf
        .join(GATK4_MUTECT2_TARGETED.out.tbi)
        .map { meta, vcf, tbi ->
            return [ meta, vcf, tbi, [] ]
        }

    GATK4_LEFTALIGNANDTRIMVARIANTS_TARGETED (
        ch_leftalignandtrimvariants_targeted,
        ch_ref_fasta.map{ meta, fasta -> [ fasta ] },
        ch_fai.map{ meta, fai -> [ fai ] },
        ch_ref_dict.map{ meta, dict -> [ dict ] }
    )
    ch_versions = ch_versions.mix(GATK4_LEFTALIGNANDTRIMVARIANTS_TARGETED.out.versions)

    ch_filtermutectcalls_targeted = GATK4_LEFTALIGNANDTRIMVARIANTS_TARGETED.out.vcf
        .join(GATK4_LEFTALIGNANDTRIMVARIANTS_TARGETED.out.tbi)
        .join(GATK4_MUTECT2_TARGETED.out.stats)
        .map { meta, vcf, tbi, stats ->
            return tuple(
                meta,
                vcf,
                tbi,
                stats,
                [],
                [],
                [],
                []
            )
        }

    GATK4_FILTERMUTECTCALLS_TARGETED (
        ch_filtermutectcalls_targeted,
        ch_ref_fasta,
        ch_fai,
        ch_ref_dict
    )
    ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS_TARGETED.out.versions)

    VARPIPE_SNPEFF_TARGETED (
        GATK4_FILTERMUTECTCALLS_TARGETED.out.vcf,
        snpeff_db,
        snpeff_cache.map { cache -> [[], cache] }
    )
    ch_versions = ch_versions.mix(VARPIPE_SNPEFF_TARGETED.out.versions)

    emit:
    ann_all           = VARPIPE_SNPEFF_ALL.out.vcf      // channel: [ val(meta), path(vcf) ]
    ann_targeted      = VARPIPE_SNPEFF_TARGETED.out.vcf // channel: [ val(meta), path(vcf) ]
    multiqc_inputs    = ch_multiqc_inputs               // channel: [ path(snpeff_report) ]
    versions          = ch_versions                     // channel: [ versions.yml ]
}
