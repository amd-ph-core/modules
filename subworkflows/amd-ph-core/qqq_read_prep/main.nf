//
// qqq's read preparation, from samplesheet runs to each sample's final reads:
//
//   1. BAM_TO_FASTQ      runs given as a BAM are converted; the layout comes from
//                        the records. Orphans are dropped, or with keep_orphans
//                        become runs of a single-end sample of their own.
//   2. FASTQ_TRIM_MERGE  Falco (optional), fastp/fastplong per run, runs merged
//   3. BBMAP_BBDUK       rRNA depletion, when an rRNA reference is given
//   4. FASTQ_DEHOST      host removal with the chosen dehoster
//   5. RASUSA_READS      downsampling: every sample when max_reads is set (rasusa
//                        writes a sample under the cap out whole), none otherwise
//
// Every intermediate the run report describes is emitted.
//
include { BAM_TO_FASTQ     } from '../bam_to_fastq/main'
include { FASTQ_TRIM_MERGE } from '../fastq_trim_merge/main'
include { FASTQ_DEHOST     } from '../fastq_dehost/main'
include { BBMAP_BBDUK      } from '../../../modules/amd-ph-core/bbmap/bbduk/main'
include { RASUSA_READS     } from '../../../modules/amd-ph-core/rasusa/reads/main'

workflow QQQ_READ_PREP {

    take:
    ch_runs             // channel: [ val(meta), val(sample), val(read_type), path(reads) ], one per samplesheet run.
                        //          reads: FASTQ(s), or one BAM (named *.bam), whose layout is read from it.
                        //          sample = [ id, run, runs ] as FASTQ_TRIM_MERGE takes it, plus for a sample with
                        //          BAM runs: bam_runs (how many of its runs are BAMs) and orphan_id (the sample
                        //          its orphans become with val_keep_orphans).
    val_keep_orphans    // value:   true to keep a paired BAM's orphans as the sample orphan_id
    val_run_falco       // value:   FASTQ_TRIM_MERGE
    val_skip_trim       // value:   FASTQ_TRIM_MERGE
    val_save_failed     // value:   FASTQ_TRIM_MERGE
    val_rrna_ref        // value:   rRNA reference FASTA for bbduk, or [] to skip rRNA depletion
    val_dehoster        // value:   FASTQ_DEHOST
    val_host_index      // value:   FASTQ_DEHOST
    val_host_index_name // value:   FASTQ_DEHOST
    val_max_reads       // value:   downsample each sample to this many reads (pairs), or 0/null for none
    val_seed            // value:   rasusa's seed

    main:
    ch_by_format = ch_runs.branch { _meta, _sample, _read_type, reads ->
        bam:   !(reads instanceof List) && reads.name.endsWith('.bam')
        fastq: true
    }

    //
    // 1. BAM input
    //
    BAM_TO_FASTQ( ch_by_format.bam.map { meta, _sample, _read_type, bam -> [ meta, bam ] } )

    // A BAM's meta gains single_end, so its sample and read type rejoin on the run id.
    ch_bam_runs = ch_by_format.bam.map { meta, sample, read_type, bam -> [ meta.id, sample, read_type, bam ] }

    ch_converted = BAM_TO_FASTQ.out.reads
        .map { meta, reads -> [ meta.id, meta, reads ] }
        .join(ch_bam_runs, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, reads, sample, read_type, _bam ->
            if (read_type == 'long' && !meta.single_end) {
                error "QQQ_READ_PREP: sample '${sample.id}' is long-read but its BAM holds read pairs. Long-read samples must be single-end."
            }
            [ meta, sample, read_type, reads ]
        }

    // Per BAM run: what the conversion found, for the run report and the orphan sample.
    ch_conversions = BAM_TO_FASTQ.out.orphans
        .map { meta, files, n_orphans -> [ meta.id, meta, files, n_orphans ] }
        .join(ch_bam_runs, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, files, n_orphans, sample, read_type, bam ->
            [ meta, sample, read_type, bam, [ layout: meta.single_end ? 'single-end' : 'paired-end', orphan_files: files, orphan_reads: n_orphans ] ]
        }

    // Kept orphans: every BAM run of a sample is collected first (bam_runs is known
    // from the samplesheet), then the runs that had any become runs of orphan_id.
    // How many that is exists only once the BAMs are read, so it is counted here.
    ch_orphan_runs = channel.empty()
    if (val_keep_orphans) {
        ch_orphan_runs = ch_conversions
            .map { _meta, sample, read_type, bam, found -> [ groupKey(sample.id, sample.bam_runs), sample, read_type, bam, found ] }
            .groupTuple(by: 0)
            .flatMap { _key, samples, read_types, bams, founds ->
                def files = samples.indices.toList()
                    .sort { a, b -> samples[a].run <=> samples[b].run }
                    .collectMany { i -> founds[i].orphan_files.collect { f -> [ f, bams[i], read_types[i] ] } }
                def orphan_id = samples[0].orphan_id
                files.withIndex().collect { fbr, i ->
                    [
                        [ id: "${orphan_id}_run${i + 1}", single_end: true ],
                        [ id: orphan_id, run: i + 1, runs: files.size() ],
                        fbr[2],
                        fbr[0],
                        fbr[1],
                    ]
                }
            }
    }

    // Every run as FASTQ (BAMs converted, kept orphans added), and what its samplesheet
    // row named: its FASTQs, or the BAM it came from.
    ch_fastq_runs = ch_by_format.fastq
        .mix(ch_converted)
        .mix(ch_orphan_runs.map { meta, sample, read_type, reads, _bam -> [ meta, sample, read_type, reads ] })
    ch_run_sources = ch_by_format.fastq.map { meta, _sample, _read_type, reads -> [ meta.id, [reads].flatten() ] }
        .mix(ch_by_format.bam.map { meta, _sample, _read_type, bam -> [ meta.id, [ bam ] ] })
        .mix(ch_orphan_runs.map { meta, _sample, _read_type, _reads, bam -> [ meta.id, [ bam ] ] })

    //
    // 2. Trim and merge
    //
    FASTQ_TRIM_MERGE(
        ch_fastq_runs,
        val_run_falco,
        val_skip_trim,
        val_save_failed
    )

    //
    // 3. rRNA depletion
    //
    ch_rrna_reads   = channel.empty()
    ch_rrna_matched = channel.empty()
    ch_rrna_stats   = channel.empty()
    ch_rrna_mqc     = channel.empty()
    ch_rrna_log     = channel.empty()
    ch_predehost    = FASTQ_TRIM_MERGE.out.reads
    if (val_rrna_ref) {
        BBMAP_BBDUK(
            FASTQ_TRIM_MERGE.out.reads.map { meta, _read_type, reads -> [ meta, reads ] },
            val_rrna_ref
        )
        ch_rrna_reads   = BBMAP_BBDUK.out.reads
        ch_rrna_matched = BBMAP_BBDUK.out.matched
        ch_rrna_stats   = BBMAP_BBDUK.out.stats
        ch_rrna_mqc     = BBMAP_BBDUK.out.stats_mqc
        ch_rrna_log     = BBMAP_BBDUK.out.log
        // bbduk drops the read type, which the dehoster needs; it rejoins on meta.
        ch_predehost = BBMAP_BBDUK.out.reads
            .join(FASTQ_TRIM_MERGE.out.reads.map { meta, read_type, _reads -> [ meta, read_type ] }, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, reads, read_type -> [ meta, read_type, reads ] }
    }

    //
    // 4. Host removal
    //
    FASTQ_DEHOST( ch_predehost, val_dehoster, val_host_index, val_host_index_name )

    //
    // 5. Downsampling: all samples or none, on the flag alone
    //
    ch_downsampled = channel.empty()
    ch_final       = FASTQ_DEHOST.out.reads
    if (val_max_reads) {
        RASUSA_READS( FASTQ_DEHOST.out.reads, val_max_reads, val_seed )
        ch_downsampled = RASUSA_READS.out.reads
        ch_final       = RASUSA_READS.out.reads
    }

    emit:
    reads          = ch_final                       // channel: [ val(meta), path(reads) ], each sample's final reads
    runs           = ch_fastq_runs                  // channel: [ val(meta), val(sample), val(read_type), path(reads) ], every run as FASTQ (BAMs converted, kept orphans added)
    run_sources    = ch_run_sources                 // channel: [ val(run_id), [ path(file) ] ], the FASTQs or BAM each run's samplesheet row named
    conversions    = ch_conversions.map { meta, sample, _read_type, bam, found -> [ meta, sample, bam, found ] }
                                                    // channel: [ val(meta), val(sample), path(bam), val([ layout, orphan_files, orphan_reads ]) ], per BAM run
    trimmed        = FASTQ_TRIM_MERGE.out.trimmed        // channel: [ val(meta), path(reads) ], per run; empty with val_skip_trim
    fastp_json     = FASTQ_TRIM_MERGE.out.fastp_json     // channel: [ val(meta), path(json) ], per short-read run
    fastplong_json = FASTQ_TRIM_MERGE.out.fastplong_json // channel: [ val(meta), path(json) ], per long-read run
    fastp_fail     = FASTQ_TRIM_MERGE.out.fastp_fail     // channel: [ val(meta), path(reads) ], per short-read run, with val_save_failed
    fastplong_fail = FASTQ_TRIM_MERGE.out.fastplong_fail // channel: [ val(meta), path(reads) ], per long-read run, with val_save_failed
    falco_txt      = FASTQ_TRIM_MERGE.out.falco_txt      // channel: [ val(meta), path(txt) ], per run, with val_run_falco
    merged         = FASTQ_TRIM_MERGE.out.merged         // channel: [ val(meta), path(reads) ], multi-run samples only
    rrna_reads     = ch_rrna_reads                  // channel: [ val(meta), path(reads) ], per sample, rRNA removed; empty without val_rrna_ref
    rrna_matched   = ch_rrna_matched                // channel: [ val(meta), path(reads) ], the reads that matched the rRNA reference
    rrna_stats     = ch_rrna_stats                  // channel: [ val(meta), path(txt) ], bbduk's match statistics
    rrna_stats_mqc = ch_rrna_mqc                    // channel: [ val(meta), path(txt) ], bbduk's summary, for MultiQC
    rrna_log       = ch_rrna_log                    // channel: [ val(meta), path(log) ], bbduk's log
    dehosted       = FASTQ_DEHOST.out.reads         // channel: [ val(meta), path(reads) ], host reads removed, or for sra_scrubber without -x kept with their bases masked as N; unchanged for none
    dehost_json    = FASTQ_DEHOST.out.json          // channel: [ val(meta), path(json) ], the dehoster's report; empty for none
    host_bam       = FASTQ_DEHOST.out.host_bam      // channel: [ val(meta), path(bam) ], hostess with --bam only
    host_bai       = FASTQ_DEHOST.out.host_bai      // channel: [ val(meta), path(bai) ], with host_bam
    host_stats     = FASTQ_DEHOST.out.host_stats    // channel: [ val(meta), path(stats) ], samtools stats of the host BAM
    host_reads     = FASTQ_DEHOST.out.host_reads    // channel: [ val(meta), path(reads) ], sra_scrubber's removed reads, with -r
    downsampled    = ch_downsampled                 // channel: [ val(meta), path(reads) ], per sample; empty without val_max_reads
}
