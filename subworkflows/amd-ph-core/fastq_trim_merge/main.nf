//
// Per-run read QC and trimming, then each sample's runs merged into one set of reads.
//
//   Falco             raw-read QC per run, when run_falco is true
//   fastp/fastplong   trimming and filtering per run (short / long reads), unless
//                     skip_trim; with save_failed the reads the filter drops are kept
//   CAT_FASTQ         a sample's runs concatenated in run order; a single-run
//                     sample passes through
//
// A run says which sample it belongs to, its place in it, and how many runs the
// sample has. The count sizes the groupKey, so each sample is merged as soon as its
// own runs are done instead of waiting on the slowest run of the whole batch.
//
include { FALCO     } from '../../../modules/amd-ph-core/falco/main'
include { FASTP     } from '../../../modules/nf-core/fastp/main'
include { FASTPLONG } from '../../../modules/amd-ph-core/fastplong/main'
include { CAT_FASTQ } from '../../../modules/nf-core/cat/fastq/main'

workflow FASTQ_TRIM_MERGE {

    take:
    ch_runs         // channel: [ val(meta), val(sample), val(read_type), path(reads) ], one per run;
                    //          sample = [ id: <merged sample id>, run: <1-based order>, runs: <run count> ],
                    //          read_type 'short' (fastp) or 'long' (fastplong)
    val_run_falco   // value:   true to run Falco on the raw reads
    val_skip_trim   // value:   true to merge the runs untrimmed
    val_save_failed // value:   true to keep the reads fastp/fastplong's filter removes

    main:
    ch_falco_txt      = channel.empty()
    ch_fastp_json     = channel.empty()
    ch_fastplong_json = channel.empty()
    ch_fastp_fail     = channel.empty()
    ch_fastplong_fail = channel.empty()
    ch_trimmed        = channel.empty()

    if (val_run_falco) {
        FALCO( ch_runs.map { meta, _sample, _read_type, reads -> [ meta, reads ] } )
        ch_falco_txt = FALCO.out.txt
    }

    ch_processed = ch_runs
    if (!val_skip_trim) {
        ch_by_type = ch_runs.branch { _meta, _sample, read_type, _reads ->
            short_reads: read_type == 'short'
            long_reads:  read_type == 'long'
            other:       true
        }
        ch_by_type.other.subscribe { meta, _sample, read_type, _reads ->
            error "FASTQ_TRIM_MERGE: run '${meta.id}' has read_type '${read_type}'; expected short or long"
        }
        FASTP(
            ch_by_type.short_reads.map { meta, _sample, _read_type, reads -> [ meta, reads, [] ] },
            false,          // discard_trimmed_pass
            val_save_failed,    // save_trimmed_fail
            false           // save_merged
        )
        FASTPLONG(
            ch_by_type.long_reads.map { meta, _sample, _read_type, reads -> [ meta, reads ] },
            [],             // adapter_fasta
            false,          // discard_trimmed_pass
            val_save_failed     // save_trimmed_fail
        )
        ch_fastp_json     = FASTP.out.json
        ch_fastplong_json = FASTPLONG.out.json
        ch_fastp_fail     = FASTP.out.reads_fail
        ch_fastplong_fail = FASTPLONG.out.reads_fail
        ch_trimmed        = FASTP.out.reads.mix(FASTPLONG.out.reads)

        // The trimmers drop the sample and read type, so they rejoin on the run's meta.
        ch_processed = ch_trimmed
            .join(ch_runs.map { meta, sample, read_type, _reads -> [ meta, sample, read_type ] }, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, reads, sample, read_type -> [ meta, sample, read_type, reads ] }
    }

    // Group each sample's runs. A sized groupKey releases a sample when its own runs
    // are done; unsized, groupTuple would hold every sample until the slowest run of
    // the batch finished.
    ch_processed
        .map { meta, sample, read_type, reads -> [ groupKey(sample.id, sample.runs), meta, sample, read_type, reads ] }
        .groupTuple(by: 0)
        .map { key, metas, samples, read_types, reads_list ->
            def sample_id = key.toString()
            if (metas*.single_end.unique().size() > 1) {
                error "FASTQ_TRIM_MERGE: sample '${sample_id}': its runs are a mix of single-end and paired-end. Multiple runs of a sample must be the same."
            }
            if (read_types.unique().size() > 1) {
                error "FASTQ_TRIM_MERGE: sample '${sample_id}': its runs are a mix of short and long reads. Multiple runs of a sample must be the same."
            }
            // In run order, not arrival order, so the merged file is reproducible.
            def ordered = samples.indices.collect { i -> [ samples[i].run, reads_list[i] ] }.sort { a, b -> a[0] <=> b[0] }
            [ [ id: sample_id, single_end: metas[0].single_end ], read_types[0], ordered.collect { run -> run[1] } ]
        }
        .branch { _meta, _read_type, runs ->
            single:   runs.size() == 1
            multiple: true
        }
        .set { ch_grouped }

    CAT_FASTQ( ch_grouped.multiple.map { meta, _read_type, runs -> [ meta, runs.flatten() ] } )

    ch_merged = CAT_FASTQ.out.reads
        .join(ch_grouped.multiple.map { meta, read_type, _runs -> [ meta, read_type ] }, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, reads, read_type -> [ meta, read_type, reads ] }
        .mix(ch_grouped.single.map { meta, read_type, runs -> [ meta, read_type, runs[0] ] })

    emit:
    reads          = ch_merged          // channel: [ val(meta), val(read_type), path(reads) ], one per sample
    merged         = CAT_FASTQ.out.reads // channel: [ val(meta), path(reads) ], multi-run samples only
    trimmed        = ch_trimmed         // channel: [ val(meta), path(reads) ], per run (empty with val_skip_trim)
    fastp_json     = ch_fastp_json      // channel: [ val(meta), path(json) ], per short-read run
    fastplong_json = ch_fastplong_json  // channel: [ val(meta), path(json) ], per long-read run
    fastp_fail     = ch_fastp_fail      // channel: [ val(meta), path(reads) ], with val_save_failed
    fastplong_fail = ch_fastplong_fail  // channel: [ val(meta), path(reads) ], with val_save_failed
    falco_txt      = ch_falco_txt       // channel: [ val(meta), path(txt) ], with val_run_falco
}
