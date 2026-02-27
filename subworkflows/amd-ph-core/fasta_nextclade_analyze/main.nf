include { NEXTCLADE_DATASETGET   } from "../../../modules/amd-ph-core/nextclade/datasetget/main"
include { NEXTCLADE_RUN          } from "../../../modules/amd-ph-core/nextclade/run/main"

workflow FASTA_NEXTCLADE_ANALYZE {

    take:
    ch_fasta          // channel: [ val(meta), path(fasta) ]
    ch_nextclade_db   // channel: [ val(dataset_name) ]
    ch_nextclade_tag  // string:  nextclade_tag

    main:

    ch_versions = Channel.empty()

    NEXTCLADE_DATASETGET(
        ch_nextclade_db,
        ch_nextclade_tag
    )
    ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)

    NEXTCLADE_RUN(
        ch_fasta,
        NEXTCLADE_DATASETGET.out.dataset
    )
    ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions)

    emit:
    csv               = NEXTCLADE_RUN.out.csv                     // channel: [ .csv ]
    tsv               = NEXTCLADE_RUN.out.tsv                     // channel: [ .tsv ]
    json              = NEXTCLADE_RUN.out.json                    // channel: [ .json ]
    json_auspice      = NEXTCLADE_RUN.out.json_auspice            // channel: [ .auspice.json ]
    ndjson            = NEXTCLADE_RUN.out.ndjson                  // channel: [ .ndjson ]
    fasta_aligned     = NEXTCLADE_RUN.out.fasta_aligned           // channel: [ .aligned.fasta ]
    fasta_translation = NEXTCLADE_RUN.out.fasta_translation       // channel: [ *_translation.*.fasta ]
    versions          = ch_versions                               // channel: [ versions.yml ]
}
