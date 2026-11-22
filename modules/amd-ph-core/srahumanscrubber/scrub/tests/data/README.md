# Test data for srahumanscrubber/scrub

- `chr22_test.w32.db`: a STAT k-mer database (k = 32, window 32; 1,245 k-mers) built with
  `srahumanscrubber/builddb` from nf-core's `genomics/homo_sapiens/genome/genome.fasta`
  (40 kb of chr22), in the same format as NCBI's `human_filter.db`.
- `human_chr22_{1,2}.fastq.gz`: 50 error-free 2 x 100 bp pairs (300 bp fragments) simulated from
  that genome with a fixed seed, so every pair is identified as human by the database above.
