# fastq_bam_mapping_scripts

First step of the sequencing analysis stage: aligns raw paired/single-end
`.fastq.gz` reads (low-coverage PGT-A embryo sequencing) to the GRCh38 reference
genome and produces one BAM file per input fastq. Everything here is written as
SLURM job scripts for the Amarel cluster.

All the actual scripts live in **[scripts/](scripts/)** — see that folder's README
for details on each one.

## Layout

- `scripts/` — bwa-mem alignment SLURM scripts (job arrays split across ~40,000
  fastq files in batches of 1000), a script to build the BAM merge job list, and a
  `test/` subfolder of earlier/scratch versions of the same scripts.
