# 4_Annotations_4_2026

# NOTE: ALL ANNOTATION WORK DONE ON LAB SERVER (minus filtered VCF creation for ANNOVAR input)

Takes the significant variants from the GLM association step (`3_GLM_4_2026`) and annotates them against genomic/functional/clinical reference databases to arrive at a list of candidate genes potentially linked to embryo aneuploidy. Everything for this stage lives under `annotation_files/` — see `annotation_files/README.md` for the full file-by-file breakdown.

## Pipeline flow

1. **Filter GLM output** (`annotation_files/glm_output/`) — the combined GLM results (`all_chromosomes_combined.csv`, from `3_GLM_4_2026`) are filtered to variants with `p < 1e-5` (111 variants), then those exact sites are pulled out of the imputed genotype VCFs (`bcftools`, on the cluster) and matched against dbSNP for authoritative rsIDs.
2. **ANNOVAR annotation** (`annotation_files/files_to_run_annovar/`) — the filtered VCF is run through ANNOVAR (`table_annovar.pl`, hg38) against RefGene/EnsGene/KnownGene gene models plus gnomAD (v2.1.1 exome, v3.1.2 genome, v4.1 exome/genome), ClinVar, dbNSFP, and dbSNP (avsnp151), producing one merged multi-database annotation table per variant.
3. **Cross-reference against external databases** (`annotation_files/`, root-level notebooks) — the annotated variants are further matched by rsID/position against:
   - the **GWAS Catalog** (prior trait associations),
   - **GTEx v10 eQTL and sQTL** data across ~50 tissues (expression/splicing effects),
   - the **snpGeneSets** R package (SNP → Entrez gene → MSigDB gene set mapping),
   - a raw **dbSNP** VCF mirror (secondary rsID verification against ANNOVAR's own avsnp151 calls).
   These results are aggregated into one table per variant (rsID-keyed), with all candidate genes pulled from every source column and converted between Ensembl and HGNC identifiers.
4. **Build the candidate gene list and per-gene table** — unique candidate genes from the aggregated variant table are extracted, optionally filtered to protein-coding only, and then enriched with **HGNC/Ensembl gene info**, **mouse ortholog mapping** (Alliance of Genome Resources + MGI), **IMPC mouse-knockout phenotype and viability data**, and **mouse oocyte expression data** (GEO GDS813), plus a manual "known meiosis gene" literature cross-reference. Several oocyte-expression cutoffs are tested; the primary output filters to genes with any oocyte expression data present.

## Key outputs

- `annotation_files/complete_filtered_aggregated_genes_df.tsv` — one row per significant variant (rsID), with every annotation source's gene calls aggregated, protein-coding genes only.
- `annotation_files/complete_variant_table_ALL_GENES_protein_and_non_protein_coding_genes.tsv` — the same, without the protein-coding filter.
- `annotation_files/candidate_ENS_gene_list.txt` / `candidate_HGNC_gene_list.txt` — the unique candidate gene list (Ensembl IDs / HGNC symbols) derived from the variant table above.
- `annotation_files/main_gene_table_HAS_EXP_filtered.csv` — **the primary candidate-gene-level output of this stage**: one row per candidate gene with Ensembl/HGNC IDs, associated variant rsIDs, mouse ortholog and IMPC viability/phenotype data, oocyte expression values, eQTL/sQTL flags, and literature cross-reference — restricted to genes with oocyte expression data. This file is also copied into `5_Gene_Set_Analysis_4_2026`, the next stage.

See `annotation_files/README.md` for the full grouped listing of every notebook and intermediate file, and each subfolder's own `README.md` for details on the reference databases used (dbSNP, GTEx, ANNOVAR's humandb, and the vendored `snpGeneSets` R package).
