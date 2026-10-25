# Data

Inputs of the analysis. `checksums.md5` lists the MD5 of each file; for a
`.gz` file it is the MD5 of the decompressed content. `make check-data`
verifies them. License: CC BY 4.0 (`LICENSE` in this folder).

Some files are *published intermediates*: tables produced by an upstream
step and used for the figure in the paper. The scripts recompute them, and
`scripts/verify.py` compares the recomputed version with the published one.

## proteomics/

| File | Content |
|---|---|
| `protein_groups.txt.gz` | MaxLFQ protein-group quantities from DIA-NN 1.8.1, 18 runs (cell extract and supernatant at 34, 37 and 40 °C, 3 replicates). Rows are protein groups; columns are the DIA files. |
| `peptides_counts.txt` | Number of distinct peptides per protein group (`Var1` protein group, `Freq` count). |
| `column_mapping.csv` | DIA file name, condition (`C34`, `C37`, `C40` cell extract; `SN34`, `SN37`, `SN40` supernatant) and replicate. |
| `all_cellular_comparisons.csv.gz` | Published intermediate. limma results for the cell extract: quantities of all 18 runs, then logFC, AveExpr, t, P.Value, adj.P.Val, B and −log10 P values for C34_C37, C34_C40 and C37_C40, two flags marking imputed proteins, and protein descriptions. |

## rnaseq/

| File | Content |
|---|---|
| `featurecounts_tb927_v68.tsv.gz` | featureCounts read-pair counts per gene against the TREU927 release 68 reference, 27 samples. The first five columns are the gene's contig, start, end, strand and length. |
| `featurecounts_bes.tsv.gz` | The same against the Tb427 bloodstream-form expression sites. |
| `GC_content_927.txt.gz`, `GC_content_BES.txt.gz` | `bedtools nuc` output for the genes of each reference; columns `4_usercol` (gene), `6_pct_gc` (GC fraction) and `13_seq_len` (length) are used. |
| `WT40_vs_WT34.txt` | Published intermediate. edgeR results, wild type 40 °C vs 34 °C: logFC, logCPM, F, PValue, FDR. |

Sample names, genotypes, temperatures and FASTQ checksums are in
`upstream/samples_rnaseq.tsv`.

## annotation/ and reference/

| File | Content |
|---|---|
| `annotation/merged_fixed_6.gtf.gz` | TREU927 annotation, chromosomes 1 to 11: TriTrypDB release 68 gene models (source `VEuPathDB`) and transcript models curated by hand in the Apollo genome browser (source `.`), which include the UTRs. |
| `reference/TriTrypDB-68_TbruceiTREU927_Genome.fasta.gz` | TriTrypDB release 68 *T. brucei* TREU927 genome. |

## sequences/

| File | Content |
|---|---|
| `trypanosoma_sequences_{5utr,cds,3utr}.csv.gz` | Columns `gene_id` and `sequence`, 6,812 genes, same row order in the three files. Recomputed exactly by `scripts/05_extract_sequences.py`. |

## embeddings/

| File | Content |
|---|---|
| `embeddings_927_hyenadna-medium-450k-seqlen-hf_trypanosoma_sequences_{5utr,cds,3utr}.csv.npy` | NumPy arrays, 6,812 × 256 float32, rows in the order of the matching sequence file. HyenaDNA medium-450k mean embeddings (`scripts/06_compute_embeddings.py`). |

## gene_universe/

| File | Content |
|---|---|
| `ml_dataset_gene_ids.txt` | The 6,923 genes of `TB/ML_DATASET.csv.gz` in github.com/mtinti/decoding-gene-expression (commit cd43e18). The Up and Down lists are restricted to these genes. |
