//
// Convert BAMs to FASTQ without knowing their layout in advance. A BAM may be
// unaligned or aligned, in any sort order; whether it is paired is read from its
// records, and comes back as meta.single_end.
//
// SAMTOOLS_COLLATEFASTQ writes read 1, read 2, other (no 0x1 flag) and singleton
// (flagged paired, mate absent) for every BAM and reports the reads in each as
// values, so the layout is decided on those counts without opening a file:
//
//   read pairs                   paired-end; other + singleton are orphans
//   no pairs, one unpaired file  single-end (usually other; singleton when every
//                                record is flagged paired but no mate is present)
//   no reads, or no pairs and both unpaired files hold reads   refused
//
// One module, on purpose. The guideline that a subworkflow chains two or more modules
// is there so a subworkflow adds something its module does not; what this one adds is
// the layout decision above, which is channel logic on the module's outputs, not
// another tool. Inlined, that logic would be copied into each pipeline that converts
// BAMs, and could not be tested on its own; as a subworkflow it is written once and
// tested here.
//
include { SAMTOOLS_COLLATEFASTQ } from '../../../modules/amd-ph-core/samtools/collatefastq/main'

workflow BAM_TO_FASTQ {

    take:
    ch_bam // channel: [ val(meta), path(bam) ]; meta.single_end is not read

    main:
    SAMTOOLS_COLLATEFASTQ (
        ch_bam,
        channel.value([ [:], [], [] ])  // no reference: BAM (or SAM) only
    )

    // Every output is emitted for every BAM with the meta it was given, so they
    // join on meta as they are.
    SAMTOOLS_COLLATEFASTQ.out.read_counts
        .join(SAMTOOLS_COLLATEFASTQ.out.fastq, failOnMismatch: true, failOnDuplicate: true)
        .join(SAMTOOLS_COLLATEFASTQ.out.fastq_other, failOnMismatch: true, failOnDuplicate: true)
        .join(SAMTOOLS_COLLATEFASTQ.out.fastq_singleton, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_bam, failOnMismatch: true, failOnDuplicate: true)
        .multiMap { meta, n_1, _n_2, n_other, n_singleton, pair, other, singleton, bam ->
            def paired   = (n_1 as long) > 0
            def unpaired = [ [ other, n_other as long ], [ singleton, n_singleton as long ] ].findAll { _f, n -> n > 0 }
            def reads
            def orphans  = []
            if (paired) {
                reads   = pair.sort { f -> f.name }
                orphans = unpaired
            }
            else if (unpaired.size() == 1) {
                reads = unpaired[0][0]
            }
            else if (unpaired.isEmpty()) {
                error "BAM_TO_FASTQ: '${meta.id}': the BAM ${bam.name} holds no primary reads."
            }
            else {
                error "BAM_TO_FASTQ: '${meta.id}': the BAM ${bam.name} holds no read pairs but both unflagged reads and reads flagged paired with no mate, so its layout is ambiguous. Convert it to FASTQ by hand."
            }
            def out_meta = meta + [ single_end: !paired ]
            reads:   [ out_meta, reads ]
            orphans: [ out_meta, orphans.collect { f, _n -> f }, orphans.sum(0L) { _f, n -> n } ]
        }
        .set { ch_converted }

    emit:
    reads   = ch_converted.reads   // channel: [ val(meta), path(reads) ], meta.single_end set from the records
    orphans = ch_converted.orphans // channel: [ val(meta), [ path(fastq) ], val(n_reads) ], one per BAM; empty list and 0 when it had none
}
