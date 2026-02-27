//
// varpipe core subworkflow
//

include { TRIMMOMATIC                                                           } from '../../../modules/amd-ph-core/trimmomatic/main'
include { VARPIPE_PROCESSBAM                                                    } from '../varpipe_processbam/main'
include { VARPIPE_COVERAGEANALYSIS                                              } from '../../../modules/amd-ph-core/varpipe/coverageanalysis/main'
include { VARPIPE_STRUCTURALVARIANTS                                            } from '../../../modules/amd-ph-core/varpipe/structuralvariants/main'
include { VARPIPE_VARIANTANALYSIS                                               } from '../varpipe_variantanalysis/main'
include { VARPIPE_CREATEANNOTATIONS as VARPIPE_CREATEANNOTATIONS_TARGETED       } from '../../../modules/amd-ph-core/varpipe/createannotations/main'
include { VARPIPE_CREATEANNOTATIONS as VARPIPE_CREATEANNOTATIONS_ALL            } from '../../../modules/amd-ph-core/varpipe/createannotations/main'
include { VARPIPE_PARSEANNOTATIONS as VARPIPE_PARSEANNOTATIONS_TARGETED         } from '../../../modules/amd-ph-core/varpipe/parseannotations/main'
include { VARPIPE_PARSEANNOTATIONS as VARPIPE_PARSEANNOTATIONS_ALL              } from '../../../modules/amd-ph-core/varpipe/parseannotations/main'
include { VARPIPE_LINEAGE                                                       } from '../../../modules/amd-ph-core/varpipe/lineage/main'
include { VARPIPE_REPORT                                                        } from '../../../modules/amd-ph-core/varpipe/report/main'
include { VARPIPE_INTERPRETATION                                                } from '../../../modules/amd-ph-core/varpipe/interpretation/main'
include { VARPIPE_PDF                                                           } from '../../../modules/amd-ph-core/varpipe/pdf/main'
include { VARPIPE_TAR                                                           } from '../../../modules/amd-ph-core/varpipe/tar/main'

workflow VARPIPE_CORE {
    take:
    ch_cleaned_reads // channel (mandatory): [ val(meta), path(reads) ]
    datestamp        // string (mandatory): timestamp for output identification
    ref_fasta        // channel (mandatory): [ path(params.ref_fasta) ]
    bedlist_amp      // channel (mandatory): [ path(params.bedlist_amp) ]
    bedlist_one      // channel (mandatory): [ path(params.bedlist_one) ]
    bedlist_two      // channel (mandatory): [ path(params.bedlist_two) ]
    bedstruct        // channel (mandatory): [ path(params.bedstruct) ]
    mutationloci     // channel (mandatory): [ path(params.mutationloci) ]
    intervals        // channel (mandatory): [ path(params.intervals) ]
    reported         // channel (mandatory): [ path(params.reported) ]
    lineages         // channel (mandatory): [ path(params.lineages) ]
    snpeff_db        // channel (mandatory): [ path(params.snpeff_db) ]
    snpeff_cache     // channel (mandatory): [ path(params.snpeff_cache) ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_inputs = Channel.empty()

    // Varpipe-specific pipeline version for stats file output and report generation
    // Maintained at subworkflow level for consistency across different workflow contexts
    def pipeline_version = '1.1.3'

    // QC TRIMMOMATIC
    TRIMMOMATIC(
        ch_cleaned_reads
    )
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)
    ch_multiqc_inputs = ch_multiqc_inputs.mix(TRIMMOMATIC.out.trimmed_reads.map { meta, reads -> reads })

    ch_ref_fasta = Channel.fromPath(ref_fasta)
        .map { fasta -> [ [ id: 'reference' ], fasta] }
        .collect()

    VARPIPE_PROCESSBAM (
        TRIMMOMATIC.out.trimmed_reads,
        ch_ref_fasta
    )
    ch_versions = ch_versions.mix(VARPIPE_PROCESSBAM.out.versions)
    ch_multiqc_inputs = ch_multiqc_inputs.mix(VARPIPE_PROCESSBAM.out.multiqc_inputs)

    ch_coverageanalysis_input = VARPIPE_PROCESSBAM.out.coverage
        .join(VARPIPE_PROCESSBAM.out.unmapped)
        .join(VARPIPE_PROCESSBAM.out.mapped)

    VARPIPE_COVERAGEANALYSIS(
        ch_coverageanalysis_input,
        bedlist_amp,
        [bedlist_one, bedlist_two],
        pipeline_version
    )
    ch_versions = ch_versions.mix(VARPIPE_COVERAGEANALYSIS.out.versions)

    VARPIPE_STRUCTURALVARIANTS(
        VARPIPE_COVERAGEANALYSIS.out.qc_passed.join(ch_coverageanalysis_input),
        VARPIPE_COVERAGEANALYSIS.out.genome_region_coverage,
        bedstruct
    )
    ch_versions = ch_versions.mix(VARPIPE_STRUCTURALVARIANTS.out.versions)

    // CALLERS
    ch_gatk_input_all = VARPIPE_COVERAGEANALYSIS.out.qc_passed
        .map { meta1, _qc -> [meta1] }
        .join(VARPIPE_PROCESSBAM.out.bam)
        .join(VARPIPE_PROCESSBAM.out.bai)
        .map { meta1, bam, bai ->
            return [meta1, bam, bai, []]
        }

    ch_gatk_input_targeted = VARPIPE_COVERAGEANALYSIS.out.qc_passed
        .map { meta1, _qc -> [meta1] }
        .join(VARPIPE_PROCESSBAM.out.bam)
        .join(VARPIPE_PROCESSBAM.out.bai)
        .map { meta1, bam, bai ->
            return [meta1, bam, bai, intervals]
        }

    VARPIPE_VARIANTANALYSIS(
        ch_gatk_input_all,
        ch_gatk_input_targeted,
        VARPIPE_PROCESSBAM.out.fai.collect(),
        VARPIPE_PROCESSBAM.out.reference_dict.collect(),
        ch_ref_fasta,
        snpeff_db,
        snpeff_cache
    )
    ch_versions = ch_versions.mix(VARPIPE_VARIANTANALYSIS.out.versions)
    ch_multiqc_inputs = ch_multiqc_inputs.mix(VARPIPE_VARIANTANALYSIS.out.multiqc_inputs)

    VARPIPE_CREATEANNOTATIONS_TARGETED(
        VARPIPE_VARIANTANALYSIS.out.ann_targeted
    )
    ch_versions = ch_versions.mix(VARPIPE_CREATEANNOTATIONS_TARGETED.out.versions.first())

    VARPIPE_CREATEANNOTATIONS_ALL(
        VARPIPE_VARIANTANALYSIS.out.ann_all
    )
    ch_versions = ch_versions.mix(VARPIPE_CREATEANNOTATIONS_ALL.out.versions.first())


    VARPIPE_PARSEANNOTATIONS_TARGETED(
        VARPIPE_CREATEANNOTATIONS_TARGETED.out.drLociAnnotation,
        mutationloci
    )
    ch_versions = ch_versions.mix(VARPIPE_PARSEANNOTATIONS_TARGETED.out.versions.first())

    VARPIPE_PARSEANNOTATIONS_ALL(
        VARPIPE_CREATEANNOTATIONS_ALL.out.fullAnnotation,
        mutationloci
    )
    ch_versions = ch_versions.mix(VARPIPE_PARSEANNOTATIONS_ALL.out.versions.first())

    // RUN VARPIPE_LINEAGE ANALYSIS
    VARPIPE_LINEAGE(
        VARPIPE_PARSEANNOTATIONS_ALL.out.fullFinalAnnotation,
        lineages
    )
    ch_versions = ch_versions.mix(VARPIPE_LINEAGE.out.versions.first())

    VARPIPE_COVERAGEANALYSIS.out.stats
        .join(VARPIPE_COVERAGEANALYSIS.out.target_region_coverage)
        .join(VARPIPE_PARSEANNOTATIONS_TARGETED.out.drLociFinalAnnotation)
        .set { ch_varpipe_report }

    // GENERATE ANALYSIS REPORT, INTERPRETATION & PDF REPORT

    VARPIPE_REPORT(
        ch_varpipe_report
    )
    ch_versions = ch_versions.mix(VARPIPE_REPORT.out.versions.first())

    VARPIPE_REPORT.out.summary
        .join(VARPIPE_STRUCTURALVARIANTS.out.structural_variants)
        .join(VARPIPE_CREATEANNOTATIONS_TARGETED.out.drLociAnnotation)
        .join(VARPIPE_COVERAGEANALYSIS.out.target_region_coverage)
        .set { ch_varpipe_interpretation }

    VARPIPE_INTERPRETATION(
        ch_varpipe_interpretation,
        reported
    )
    ch_versions = ch_versions.mix(VARPIPE_INTERPRETATION.out.versions.first())

    VARPIPE_PDF(
        VARPIPE_INTERPRETATION.out.summary
    )
    ch_versions = ch_versions.mix(VARPIPE_PDF.out.versions.first())

    ch_varpipe_outputs = Channel.empty()
        .mix(
            VARPIPE_COVERAGEANALYSIS.out.stats,
            VARPIPE_COVERAGEANALYSIS.out.qc_failed,
            VARPIPE_COVERAGEANALYSIS.out.target_region_coverage,
            VARPIPE_COVERAGEANALYSIS.out.genome_region_coverage,
            VARPIPE_STRUCTURALVARIANTS.out.structural_variants,
            VARPIPE_LINEAGE.out.lineage,
            VARPIPE_LINEAGE.out.lineage_report,
            VARPIPE_INTERPRETATION.out.summary,
            VARPIPE_INTERPRETATION.out.interpretation,
            VARPIPE_PDF.out.report,
            VARPIPE_PARSEANNOTATIONS_TARGETED.out.drLociFinalAnnotation,
            VARPIPE_CREATEANNOTATIONS_TARGETED.out.drLociAnnotation,
            VARPIPE_VARIANTANALYSIS.out.ann_targeted,
            VARPIPE_VARIANTANALYSIS.out.ann_all,
            VARPIPE_PARSEANNOTATIONS_ALL.out.fullFinalAnnotation,
            VARPIPE_CREATEANNOTATIONS_ALL.out.fullAnnotation
        )
        .groupTuple()
        .multiMap { meta, files ->
            meta: [ id: "${meta.id}", filenames: files.collect { it.getName() } ]
            files: files
        }

    VARPIPE_TAR(
        ch_varpipe_outputs.meta.collect(),
        ch_varpipe_outputs.files.collect(),
        datestamp
    )
    ch_versions = ch_versions.mix(VARPIPE_TAR.out.versions)

    emit:
    lineage                         = VARPIPE_LINEAGE.out.lineage_report            // channel (mandatory): [ val(meta), path(lineage_report) ]
    multiqc_inputs                  = ch_multiqc_inputs                             // channel (mandatory): [ varpipe multiqc files ]
    versions                        = ch_versions                                   // channel (mandatory): [ versions.yml ]
}
