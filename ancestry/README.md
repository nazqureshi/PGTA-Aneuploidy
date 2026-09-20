# scripts (Ancestry_inference_Scripts/scripts)

LASER-based ancestry inference pipeline. All scripts are SLURM jobs for Amarel and
use paths relative to this directory (`../data/...`, `../logs/...`), which are not
present in this checkout — they exist on the cluster at runtime. Steps are lettered
A–G in the scripts' own comments and should be run in that order.

- **`prepare_reference.sh`** (Step A/B) — Filters the 1000 Genomes phased VCFs
  (`ALL.chr<N>.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz`) down
  to biallelic SNPs with minor allele frequency > 10% per chromosome
  (`bcftools view`, run 8-way parallel across chromosomes), then concatenates all 22
  filtered chromosome VCFs into one reference VCF (`1kg_maf10.vcf.gz`) and removes
  the per-chromosome intermediates.

- **`convert_reference.sh`** (Step C/D) — Converts the merged reference VCF to
  LASER's native `.geno` format using LASER's own `vcf2geno` tool, and derives a
  `.bed` file of SNP positions from the resulting `.site` file (for use in generating
  pileups at the same sites).

- **`generate_pileups.sh`** (Step E) — For each patient BAM (listed in
  `../data/bam/bam_list.txt`), runs `samtools mpileup` restricted to the reference
  SNP positions (`1kg.bed`) against the hg38 reference genome, producing one
  `.pileup` file per sample (8-way parallel), then writes a list of all pileup files.

- **`pileupsToSeq.sh`** — Converts the per-sample pileups into LASER's `.seq` format
  using LASER's bundled `pileup2seq/pileup2seq.py` (Python 2), referencing the same
  reference `.site` file, producing `target_samples.seq`.

- **`PCA_step.sh`** (Step F) — Runs LASER's `laser` binary in PCA mode (`-pca 1`) on
  the reference genotypes (`1kg.geno`) to compute the reference PCA coordinates
  (`laser_PCA.RefPC.coord`), 8 PCs.

- **`infer_ancestry.sh`** (Step G) — Runs `laser` in its default ancestry-estimation
  mode, projecting each patient's `.seq` data onto the reference PCA
  (`-c laser_PCA.RefPC.coord -g 1kg.geno -s target_samples.seq`), 8 PCs, 5 repeated
  runs per sample (`-r 5`), 24 threads. This is the final step producing per-patient
  ancestry PC coordinates.

- **`pca.py`** — Standalone helper (not a SLURM job) that merges LASER's reference
  PCA output (`laser.RefPC.coord`) with a 1000 Genomes phenotype/population table to
  add a `popID` (superpopulation code) column, for downstream plotting/labeling.
  Note: hardcodes paths under an older/different location than the rest of this
  pipeline, so may need updated paths before reuse.

- **`laser.conf`** — LASER's parameter file template (all fields commented out /
  unset) documenting every option the `laser` binary accepts; the pipeline scripts
  above pass parameters via command-line flags instead of using this file directly.

- **`laser.log`**, **`laser.SeqPC.coord`**, **`laser.SeqPC.coord.sd`** — Output
  artifacts from a completed `infer_ancestry.sh` run (5616 individuals, 2548
  reference individuals, 8 PCs): the run log, per-sample projected PC coordinates,
  and their standard deviations across the 5 repeated runs.

## Vendored tool

- **`LASER-2.04/`** — Vendored copy of the LASER ancestry-inference tool (v2.04), see
  `LASER-2.04/ReadMe.txt` for upstream docs. Not documented further here; scripts
  above reference its `laser` and `vcf2geno` binaries and `pileup2seq/pileup2seq.py`
  script directly.
