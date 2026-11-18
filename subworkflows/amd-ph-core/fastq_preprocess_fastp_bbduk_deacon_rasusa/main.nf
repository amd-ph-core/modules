//
// Read preprocessing as QQQ runs it: trim, deplete rRNA, remove host, and optionally cap depth.
//
//   1. fastp / fastplong  adapter trimming, and every read-level QC metric      always
//   2. bbduk              rRNA depletion against a reference panel              on by default
//   3. deacon             host removal against a minimizer index               on by default
//   4. rasusa             downsampling to a fixed number of reads or pairs      off by default
//
// The ORDER is the point of this subworkflow, and both choices in it are measured.
//
// rRNA depletion runs BEFORE host removal. rRNA is the most conserved sequence there is, so a
// host index removes more than host rRNA: deacon drops a read on shared minimizers, and on
// simulated reads it removed 51-53% of fungal 18S (Candida albicans, Aspergillus fumigatus), 10%
// of Plasmodium falciparum 18S and 9% of nematode 18S, keeping none of it. Run first, bbduk takes
// all of that rRNA, from every domain, into its matched stream, where it can still be
// characterised. It also makes the two removal figures mean what they say: "rRNA" is rRNA of any
// origin, human included, and "host" is host reads that are not rRNA. The cost is that bbduk
// also reads the host reads: on a 198-sample arbovirus mNGS run, an estimated 0.38 CPU-hours on
// top of 2.98, at most 2 minutes for a 13.6 Gbase library that was 79% host.
//
// Downsampling runs LAST, after both removals. Capped earlier, the budget is spent on reads that
// are about to be discarded: on an RNA library rRNA is most of the sample -- 94% of the
// wastewater library this was built against -- so a 1M-read cap ahead of depletion leaves ~60k
// reads of anything else.
//
// The rRNA-matched stream is emitted rather than dropped. It is the input to rRNA
// characterisation, and the only place an over-aggressive panel shows up as reads rather than as
// a quietly smaller library. Because it is taken before host removal it CONTAINS HOST SEQUENCE:
// host rRNA, mitochondrial 12S/16S, and the mate of any read that matched, since bbduk moves the
// pair. Treat it as human-derived data.
//
// There is no versions channel. Every module here emits its version to the `versions` topic,
// which the pipeline collects with channel.topic('versions'); a subworkflow-level channel would
// report each tool twice.
//

include { FASTP         } from '../../../modules/nf-core/fastp/main'
include { FASTPLONG     } from '../../../modules/nf-core/fastplong/main'
include { BBMAP_BBDUK   } from '../../../modules/amd-ph-core/bbmap/bbduk/main'
include { DEACON_FILTER } from '../../../modules/amd-ph-core/deacon/filter/main'
include { RASUSA_READS  } from '../../../modules/amd-ph-core/rasusa/reads/main'

workflow FASTQ_PREPROCESS_FASTP_BBDUK_DEACON_RASUSA {

    take:
    ch_short_reads          // channel: [ val(meta), [ path(reads) ] ]  short reads, single- or paired-end
    ch_long_reads           // channel: [ val(meta), path(reads) ]      long reads, one file per sample
    ch_rrna_ref             // channel: path(fasta)  rRNA reference     -- ignored when skip_rrna_depletion
    ch_host_index           // channel: path(idx)    deacon index       -- ignored when skip_dehost
    skip_rrna_depletion     // boolean: skip bbduk
    skip_dehost             // boolean: skip deacon
    skip_downsample         // boolean: skip rasusa
    max_reads               // integer: reads (single-end) or read pairs (paired-end) to keep per sample
    downsample_seed         // integer: rasusa random seed

    main:
    // The switches are plain booleans, not channels: `if (channel)` is always true, so a step
    // switched off through a channel would silently run anyway.
    [ skip_rrna_depletion: skip_rrna_depletion, skip_dehost: skip_dehost, skip_downsample: skip_downsample ]
        .each { name, value ->
            if (!(value instanceof Boolean)) {
                error("FASTQ_PREPROCESS_FASTP_BBDUK_DEACON_RASUSA: ${name} must be a boolean, got ${value?.getClass()?.simpleName}: ${value}")
            }
        }
    if (!skip_downsample && !(max_reads instanceof Number && max_reads > 0)) {
        error("FASTQ_PREPROCESS_FASTP_BBDUK_DEACON_RASUSA: downsampling is enabled, so max_reads must be a positive integer, got: ${max_reads}")
    }

    ch_multiqc_files = channel.empty()

    //
    // STEP 1: adapter trimming. fastp's JSON report is the source of every read-level QC metric
    // -- counts, lengths, base quality, GC, duplication, adapter content -- so no separate QC
    // tool runs. What is trimmed or filtered is set entirely through ext.args by the pipeline.
    //
    FASTP(
        ch_short_reads.map { meta, reads -> [ meta, reads, [] ] },
        false,  // discard_trimmed_pass
        false,  // save_trimmed_fail
        false,  // save_merged
    )
    FASTPLONG(
        ch_long_reads,
        [],     // adapter_fasta
        false,  // discard_trimmed_pass
        false,  // save_trimmed_fail
    )
    ch_trimmed_reads = FASTP.out.reads.mix(FASTPLONG.out.reads)
    ch_multiqc_files = ch_multiqc_files.mix(
        FASTP.out.json.map { _meta, json -> json },
        FASTPLONG.out.json.map { _meta, json -> json },
    )

    //
    // STEP 2: rRNA depletion.
    //
    // `.first()` turns the reference into a value channel, so that it is broadcast to every
    // sample instead of being consumed by the first one.
    //
    ch_rrna_reads  = channel.empty()
    ch_bbduk_stats = channel.empty()
    ch_bbduk_log   = channel.empty()
    if (skip_rrna_depletion) {
        ch_depleted_reads = ch_trimmed_reads
    }
    else {
        BBMAP_BBDUK(ch_trimmed_reads, ch_rrna_ref.first())
        ch_depleted_reads = BBMAP_BBDUK.out.reads
        ch_rrna_reads     = BBMAP_BBDUK.out.matched
        ch_bbduk_stats    = BBMAP_BBDUK.out.stats
        // The exact count of records removed. #Matched in the stats file undercounts paired-end
        // removal, because one matching mate moves the pair.
        ch_bbduk_log      = BBMAP_BBDUK.out.log
        // The truncated stats copy, not the full one: MultiQC writes every row it parses into
        // its own outputs, and against SILVA SSU+LSU NR99 the full file has 93,601 rows, 23 MB
        // of multiqc_data.json per sample. The truncated copy gives the same summary from 38 KB.
        ch_multiqc_files  = ch_multiqc_files.mix(BBMAP_BBDUK.out.stats_mqc.map { _meta, stats -> stats })
    }

    //
    // STEP 3: host removal.
    //
    ch_deacon_json = channel.empty()
    if (skip_dehost) {
        ch_dehosted_reads = ch_depleted_reads
    }
    else {
        DEACON_FILTER(ch_depleted_reads, ch_host_index.first())
        ch_dehosted_reads = DEACON_FILTER.out.fastq
        ch_deacon_json    = DEACON_FILTER.out.json
        ch_multiqc_files  = ch_multiqc_files.mix(DEACON_FILTER.out.json.map { _meta, json -> json })
    }

    //
    // STEP 4: downsampling. Every sample goes through rasusa, which passes a sample already
    // under the cap through whole, so the result does not depend on a read count measured
    // upstream. That costs one pass over the reads, and only when the step is switched on.
    //
    if (skip_downsample) {
        ch_clean_reads = ch_dehosted_reads
    }
    else {
        RASUSA_READS(ch_dehosted_reads, max_reads, downsample_seed)
        ch_clean_reads = RASUSA_READS.out.reads
    }

    emit:
    reads          = ch_clean_reads          // channel: [ val(meta), [ path(reads) ] ]  trimmed, depleted, dehosted, capped
    rrna_reads     = ch_rrna_reads           // channel: [ val(meta), [ path(reads) ] ]  bbduk matched stream; contains host
    fastp_json     = FASTP.out.json          // channel: [ val(meta), path(json) ]
    fastp_html     = FASTP.out.html          // channel: [ val(meta), path(html) ]
    fastplong_json = FASTPLONG.out.json      // channel: [ val(meta), path(json) ]
    fastplong_html = FASTPLONG.out.html      // channel: [ val(meta), path(html) ]
    bbduk_stats    = ch_bbduk_stats          // channel: [ val(meta), path(txt) ]
    bbduk_log      = ch_bbduk_log            // channel: [ val(meta), path(log) ]
    deacon_json    = ch_deacon_json          // channel: [ val(meta), path(json) ]
    multiqc_files  = ch_multiqc_files        // channel: path(file)
}
