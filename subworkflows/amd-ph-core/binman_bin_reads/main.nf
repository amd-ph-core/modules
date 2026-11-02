//
// Sort reads into family-level bins: kraken2 rank-direct binning, then recruit the residual.
//
// The cascade shape matters and is measured. Each tier runs ONLY on what the previous tiers
// left, and cumulative recruitment equals the union of all tiers exactly -- residual-only
// scanning never discards a read, it just stops re-examining ones already placed. On a
// 200,000-pair subsample of an rRNA-depleted wastewater library:
//
//   stage                   scans     recruits   yield   residual left
//   kraken2 family        200,000       45,786   22.9%         154,214
//   + mmseqs translated   154,214       71,949   67.5%          82,265
//   + dense panel          82,265          841    1.0%          81,424
//
// (the mmseqs and dense-panel rows are reordered from the measured cascade to reflect the
// recommended order, which runs the cheaper, higher-yield protein tier first)
//
// Assembly-based recruitment is deliberately NOT a tier here. PenguiN recruited the most reads
// of anything tested (88.2% of the pool, 16,434 of them unique), but its unique reads are
// circular evidence: the contig was assembled from this very pool, so a read mapping back to a
// contig it helped build is not independent confirmation the way a reference hit is. Assembly
// belongs downstream of binning, operating on bins.
//
// Tier 3 is MINIMAP2_ALIGN rather than the rammap module this subworkflow carried in the binman
// pipeline. rammap is a pure-Rust reimplementation of minimap2 that we build ourselves, and for
// this tier the two are not merely similar: on 300,000 reads of a real depleted wastewater
// library against a 3,138-record dense panel spanning Picornaviridae, Caliciviridae and
// Virgaviridae, `-f 0` at both `-x lr` and `-x sr`, the two tools produced BYTE-IDENTICAL PAF
// (md5 a142379855c89399fd1e1dc9540b21ef, 314,760 records) and therefore identical read->family
// assignments (30,395 reads recruited). rammap's ouroboros patches -- `--filter-chimera` and
// `local_sw` -- are what minimap2 cannot provide, and this tier uses neither. Using the nf-core
// module instead buys a public, maintained container in place of an ECR image, which is what
// lets this subworkflow be offered to pipelines outside our account.
//
// One behaviour difference to know about: MINIMAP2_ALIGN writes an UNCOMPRESSED `.paf`, whereas
// the rammap module streamed into gzip because the combined PAF measured 974 MB for 11 families
// on one library and would plausibly reach several GB at 214. Nothing downstream cares --
// BINMAN_ASSIGNBINS reads plain or gzipped -- but the work directory does.
//
// Tier 2 is likewise the nf-core module rather than the MMSEQS_TRANSLATEDSEARCH module this
// subworkflow carried at first, and in this case the public module is not merely equivalent but
// CORRECTER. Measured on 2,000 real residual reads against the 2,622,400-ORF production target
// database: both produce the same 65,108 rows with identical query/target pairs, alignment
// coordinates, alignment lengths, e-values and bit scores. The only column that differed was
// percent identity, and ours was the wrong one -- `mmseqs easy-search` passes
// `--alignment-mode 3` to its internal `search`, which computes real sequence identity, whereas
// a bare `mmseqs search` defaults to alignment mode 2 (score and end positions only) and
// `convertalis` then reports an estimate. Re-running the search with `-a 1` to get an exact
// backtrace settles it: 22 identities over 34 aligned residues, 64.700%, exactly what
// easy-search reports, where the hand-rolled chain said 70.100%. Adding `--alignment-mode 3`
// to the bespoke module made its output byte-identical to this one, which is both the proof of
// cause and the reason not to keep maintaining a second copy.
//
// Nothing downstream read that column -- BINMAN_ASSIGNBINS uses only the query and target names
// -- so no measurement is invalidated. Anyone who filtered tier-2 hits on identity, or on
// `--min-seq-id`, was filtering on the estimate.
//
// MMSEQS_EASYSEARCH takes the residual mates as they are: `mmseqs createdb` reads gzipped FASTQ
// directly, both mates in one query database, so no FASTQ-to-FASTA step is needed. It writes an
// uncompressed `.tsv` in the 12-column BLAST format by default; set `--format-output` through
// `ext.args` if you want the narrower table. BINMAN_ASSIGNBINS reads either.
//

include { KRAKEN2HUNTER_CLASSIFY } from '../../../modules/amd-ph-core/kraken2hunter/classify/main'
include { BINMAN_RESIDUAL        } from '../../../modules/amd-ph-core/binman/residual/main'
include { MMSEQS_EASYSEARCH      } from '../../../modules/nf-core/mmseqs/easysearch/main'
include { MINIMAP2_ALIGN         } from '../../../modules/nf-core/minimap2/align/main'
include { BINMAN_ASSIGNBINS      } from '../../../modules/amd-ph-core/binman/assignbins/main'

workflow BINMAN_BIN_READS {

    take:
    ch_reads          // channel: [ val(meta), [ reads ] ]
    ch_kraken2_db     // channel: path(db)
    bin_rank          // val: rank bins resolve to, e.g. 'family'
    ch_mmseqs_db      // channel: path(db)            -- ignored unless run_mmseqs
    ch_panel          // channel: path(fasta)         -- ignored unless run_minimap2
    ch_panel_blocks   // channel: path(tsv)           -- ignored unless run_minimap2
    run_mmseqs        // val: enable the translated-search tier
    run_minimap2      // val: enable the dense-panel nucleotide tier

    main:
    // No ch_versions plumbing: every module in this subworkflow emits its version to the
    // `versions` topic, which the pipeline collects globally via channel.topic("versions").
    // Threading a versions channel through here as well would double-report each tool.
    //
    // The two tiers are switched by plain `val` booleans rather than read from `params`, which
    // is what promotion to a shared repo requires -- a consumer that has never declared
    // `params.enrich_mmseqs` would fail at evaluation. It also keeps the guard explicit: in the
    // pipeline this came from, enabling a tier by flag alone while leaving its reference
    // channel empty did not fail. It declared the process with nothing to run, so the tier was
    // skipped in silence. A caller must decide, and should decide on BOTH its flag and whether
    // it actually has the reference.

    //
    // TIER 1: classify and bin at the requested rank in one pass
    //
    KRAKEN2HUNTER_CLASSIFY(ch_reads, ch_kraken2_db, bin_rank)

    //
    // Assemble the residual pool. Both inputs are optional on the tool side, so they are
    // defaulted to an empty list per sample rather than joined on presence -- a sample whose
    // unresolved bin is empty must still produce a residual pool, because its unclassified
    // stream is typically 60-80% of the library.
    //
    ch_residual_in = ch_reads
        .map { meta, _reads -> [ meta, [] ] }
        .join(KRAKEN2HUNTER_CLASSIFY.out.unresolved,   remainder: true)
        .join(KRAKEN2HUNTER_CLASSIFY.out.unclassified, remainder: true)
        .map { meta, _placeholder, unres, unc -> [ meta, unres ?: [], unc ?: [] ] }

    BINMAN_RESIDUAL(ch_residual_in)

    //
    // TIER 2: translated search. Highest yield per unit compute of anything measured.
    //
    // When a tier is off, its channel is a per-sample PLACEHOLDER rather than Channel.empty().
    // This is the fix for a real bug: joining against an empty channel, even with
    // `remainder: true`, dropped samples -- only 1 of 3 reached BINMAN_ASSIGNBINS. A
    // placeholder keyed by meta means the downstream join is a plain inner join that always
    // has a partner, which behaves identically whether or not the tier ran.
    //
    ch_mmseqs_hits = ch_reads.map { meta, _reads -> [ meta, [] ] }
    if (run_mmseqs) {
        // The database arrives as a bare path so that this subworkflow's input signature stays
        // about the reference rather than about the module consuming it; MMSEQS_EASYSEARCH wants
        // a meta on that channel, so one is attached here. `.first()` makes it a value channel,
        // for the same reason it is needed on ch_panel_blocks below: a queue channel is consumed
        // rather than broadcast and would starve every sample after the first.
        MMSEQS_EASYSEARCH(
            BINMAN_RESIDUAL.out.residual,
            ch_mmseqs_db.map { db -> [ [ id:'target' ], db ] }.first(),
        )
        ch_mmseqs_hits = MMSEQS_EASYSEARCH.out.tsv
    }

    //
    // TIER 3: dense-panel nucleotide recruitment.
    //
    // The caller sets the preset and `-f 0` through `ext.args`. `-f 0` is not optional on a
    // dense panel: the mid-occ minimizer cap discards exactly the informative minimizers, and
    // it is what made this tier recruit 24.5% of a library where blastn found 68.4% at the same
    // threshold. There is deliberately no MAPQ filter anywhere downstream either -- on a dense
    // panel the competing alignments are other variants of the same organism, so MAPQ is near
    // zero for precisely the reads the panel exists to catch, and `-q 20` deleted 98.5% of
    // on-target reads in the arbovirus work.
    //
    ch_paf = ch_reads.map { meta, _reads -> [ meta, [] ] }
    if (run_minimap2) {
        MINIMAP2_ALIGN(
            BINMAN_RESIDUAL.out.residual,
            ch_panel.map { panel -> [ [ id:'panel' ], panel ] }.first(),
            false,  // bam_format
            '',     // bam_index_extension
            false,  // cigar_paf_format -- the CIGAR is not read; only match/block counts are
            false,  // cigar_bam
        )
        ch_paf = MINIMAP2_ALIGN.out.paf
    }

    //
    // Fold enrichment hits into the tier-1 bins.
    //
    // Plain inner joins, because every optional input above is guaranteed exactly one entry per
    // sample -- a disabled tier contributes a per-sample [] placeholder rather than
    // Channel.empty(). That is deliberate: `join(remainder: true)` against a genuinely empty
    // channel dropped samples, leaving only 1 of 3 here, and chaining two such joins made the
    // tuple arity inconsistent on top of that. Keeping both sides one-entry-per-sample means
    // the join behaves identically whether a tier ran or not, with no remainder semantics to
    // reason about.
    //
    // Written out rather than factored into a helper closure because Nextflow 26's strict
    // parser does not allow calling a locally-assigned closure as a function in a workflow body.
    ch_assign = BINMAN_RESIDUAL.out.residual
        .join(KRAKEN2HUNTER_CLASSIFY.out.bin_manifest)
        .join(ch_mmseqs_hits)
        .join(ch_paf)

    // `.first()` converts the blocks reference to a VALUE channel, and that is load-bearing.
    // A queue channel is consumed rather than broadcast, so a single-item queue channel --
    // exactly what `channel.empty().ifEmpty([])` produces -- satisfies the first sample and
    // then starves the rest. The process ran once for the whole run, publishing enriched bins
    // for one sample and silently nothing for the other two.
    BINMAN_ASSIGNBINS(ch_assign, ch_panel_blocks.ifEmpty([]).first())

    emit:
    report        = KRAKEN2HUNTER_CLASSIFY.out.report       // channel: [ meta, txt ]
    bins          = KRAKEN2HUNTER_CLASSIFY.out.bins         // channel: [ meta, [ fastq.gz ] ]
    bin_manifest  = KRAKEN2HUNTER_CLASSIFY.out.bin_manifest // channel: [ meta, tsv ]
    bins_enriched = BINMAN_ASSIGNBINS.out.bins              // channel: [ meta, [ fastq.gz ] ]
    residual      = BINMAN_ASSIGNBINS.out.residual          // channel: [ meta, [ fastq.gz ] ]
    assignment    = BINMAN_ASSIGNBINS.out.report            // channel: [ meta, tsv ]
}
