//
// Remove host reads with one of four dehosters, building the host index from a
// FASTA when one is given.
//
//   deacon        minimizer matching against a deacon index
//   sra_scrubber  NCBI's SRA Human Read Scrubber: STAT k-mers, reads removed
//                 (or masked, through SRAHUMANSCRUBBER_SCRUB's ext.args)
//   hostess       alignment, keeping the removed reads as a BAM if asked
//   hostile       alignment
//   none          reads pass through unchanged
//
// There is one process per dehoster, whatever the input's size: resources that
// should grow with the input belong in the pipeline's config, as closures over
// the task's own reads, not in a second copy of the process. The aligning
// dehosters run as _SHORT (bowtie2) and _LONG (minimap2), which differ in the
// aligner and the index they are handed, not just in resources.
//
include { DEACON_INDEX                                 } from '../../../modules/amd-ph-core/deacon/index/main'
include { DEACON_FILTER                                } from '../../../modules/amd-ph-core/deacon/filter/main'
include { SRAHUMANSCRUBBER_BUILDDB                     } from '../../../modules/amd-ph-core/srahumanscrubber/builddb/main'
include { SRAHUMANSCRUBBER_SCRUB                       } from '../../../modules/amd-ph-core/srahumanscrubber/scrub/main'
include { BOWTIE2_BUILD                                } from '../../../modules/nf-core/bowtie2/build/main'
include { HOSTESS as HOSTESS_SHORT                     } from '../../../modules/amd-ph-core/hostess/main'
include { HOSTESS as HOSTESS_LONG                      } from '../../../modules/amd-ph-core/hostess/main'
include { HOSTILE_CLEAN as HOSTILE_CLEAN_SHORT         } from '../../../modules/amd-ph-core/hostile/clean/main'
include { HOSTILE_CLEAN as HOSTILE_CLEAN_LONG          } from '../../../modules/amd-ph-core/hostile/clean/main'

workflow FASTQ_DEHOST {

    take:
    ch_reads        // channel: [ val(meta), val(read_type), path(reads) ], read_type 'short' or 'long'
    dehoster        // value:   deacon | sra_scrubber | hostess | hostile | none
    host_fasta      // value:   a host FASTA to build the index from, or [] to use host_index
    host_index      // value:   a prebuilt index for the dehoster: deacon .idx, scrubber .db, or the
                    //          hostile/hostess index directory; [] when host_fasta is given
    host_index_name // value:   hostile/hostess: the index name inside host_index (bowtie2 prefix); '' otherwise
    scrubber_window // value:   sra_scrubber: the window SRAHUMANSCRUBBER_BUILDDB builds with

    main:
    if (!(dehoster in ['deacon', 'sra_scrubber', 'hostess', 'hostile', 'none'])) {
        error "FASTQ_DEHOST: unknown dehoster '${dehoster}'; expected deacon, sra_scrubber, hostess, hostile or none"
    }
    if (dehoster != 'none' && !host_fasta && !host_index) {
        error "FASTQ_DEHOST: dehoster '${dehoster}' needs either a host FASTA to build an index from or a prebuilt index"
    }

    ch_dehosted   = channel.empty()
    ch_json       = channel.empty()
    ch_host_bam   = channel.empty()
    ch_host_bai   = channel.empty()
    ch_host_stats = channel.empty()
    ch_host_reads = channel.empty()
    ch_fastq      = ch_reads.map { meta, _read_type, reads -> [ meta, reads ] }

    if (dehoster == 'none') {
        ch_dehosted = ch_fastq
    }
    else if (dehoster == 'deacon') {
        if (host_fasta) {
            DEACON_INDEX( channel.value([ [ id: file(host_fasta).baseName ], host_fasta ]) )
            ch_index = DEACON_INDEX.out.index.map { _meta, index -> index }.first()
        }
        else {
            ch_index = channel.value(host_index)
        }
        DEACON_FILTER( ch_fastq, ch_index )
        ch_dehosted = DEACON_FILTER.out.fastq
        ch_json     = DEACON_FILTER.out.json
    }
    else if (dehoster == 'sra_scrubber') {
        if (host_fasta) {
            SRAHUMANSCRUBBER_BUILDDB( channel.value([ [ id: file(host_fasta).baseName ], host_fasta ]), [], scrubber_window )
            ch_db = SRAHUMANSCRUBBER_BUILDDB.out.db.map { _meta, db -> db }.first()
        }
        else {
            ch_db = channel.value(host_index)
        }
        SRAHUMANSCRUBBER_SCRUB( ch_fastq, ch_db )
        ch_dehosted   = SRAHUMANSCRUBBER_SCRUB.out.fastq
        ch_json       = SRAHUMANSCRUBBER_SCRUB.out.json
        ch_host_reads = SRAHUMANSCRUBBER_SCRUB.out.removed
    }
    else {
        ch_by_type = ch_reads.branch { _meta, read_type, _reads ->
            short_reads: read_type == 'short'
            long_reads:  read_type == 'long'
            other:       true
        }
        ch_by_type.other.subscribe { meta, read_type, _reads ->
            error "FASTQ_DEHOST: sample '${meta.id}' has read_type '${read_type}'; expected short or long"
        }
        ch_short = ch_by_type.short_reads.map { meta, _read_type, reads -> [ meta, reads ] }
        ch_long  = ch_by_type.long_reads.map  { meta, _read_type, reads -> [ meta, reads ] }

        // bowtie2 takes an index name inside a directory; minimap2 takes the FASTA
        // (or a directory holding <name>.mmi), which is why the two references differ.
        // The bowtie2 index is built only when some sample is short-read: built from
        // the first short-read sample's arrival, so a long-read-only run (or an empty
        // one) never builds a human-genome index it would not use.
        if (host_fasta) {
            BOWTIE2_BUILD( ch_short.first().map { _meta, _reads -> [ [ id: file(host_fasta).baseName ], host_fasta ] } )
            ch_ref_short = BOWTIE2_BUILD.out.index.map { meta, index -> [ meta.id, index ] }.first()
            ch_ref_long  = channel.value([ file(host_fasta).name, host_fasta ])
        }
        else {
            ch_ref_short = channel.value([ host_index_name, host_index ])
            ch_ref_long  = channel.value([ host_index_name, host_index ])
        }

        if (dehoster == 'hostess') {
            HOSTESS_SHORT( ch_short, ch_ref_short, 'bowtie2' )
            HOSTESS_LONG ( ch_long,  ch_ref_long,  'minimap2' )
            ch_dehosted   = HOSTESS_SHORT.out.fastq.mix(HOSTESS_LONG.out.fastq)
            ch_json       = HOSTESS_SHORT.out.json.mix(HOSTESS_LONG.out.json)
            ch_host_bam   = HOSTESS_SHORT.out.bam.mix(HOSTESS_LONG.out.bam)
            ch_host_bai   = HOSTESS_SHORT.out.bai.mix(HOSTESS_LONG.out.bai)
            ch_host_stats = HOSTESS_SHORT.out.stats.mix(HOSTESS_LONG.out.stats)
        }
        else {
            HOSTILE_CLEAN_SHORT( ch_short, ch_ref_short, 'bowtie2' )
            HOSTILE_CLEAN_LONG ( ch_long,  ch_ref_long,  'minimap2' )
            ch_dehosted = HOSTILE_CLEAN_SHORT.out.fastq.mix(HOSTILE_CLEAN_LONG.out.fastq)
            ch_json     = HOSTILE_CLEAN_SHORT.out.json.mix(HOSTILE_CLEAN_LONG.out.json)
        }
    }

    emit:
    reads      = ch_dehosted   // channel: [ val(meta), path(reads) ], host reads removed (unchanged for none)
    json       = ch_json       // channel: [ val(meta), path(json) ], the dehoster's report; empty for none
    host_bam   = ch_host_bam   // channel: [ val(meta), path(bam) ], hostess with --bam only
    host_bai   = ch_host_bai   // channel: [ val(meta), path(bai) ]
    host_stats = ch_host_stats // channel: [ val(meta), path(stats) ], hostess samtools stats of the host BAM
    host_reads = ch_host_reads // channel: [ val(meta), path(reads) ], sra_scrubber's removed reads, with -r
}
