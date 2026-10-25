# Methods

This document describes each step of the analysis with its parameters, in
the order the data flow. The script for each step is in `scripts/` or
`upstream/`.

## Samples

Bloodstream-form *T. brucei* Lister 427 cells, wild type and two independent
*ago1* null lines (A and D), were grown at 34, 37 and 40 °C for 6 hours, in
biological triplicate. RNA-seq covers
all 27 samples (NCBI SRA PRJNA1381115; runs in `upstream/samples_rnaseq.tsv`). Proteomics covers wild-type
cells at the three temperatures, as cell extract and culture supernatant, in
triplicate (PRIDE PXD071951). The figure compares 40 °C with 34 °C: the cell extract for
proteomics, and the wild-type samples for RNA-seq.

## Proteomics

### Search and quantification (`upstream/`)

DIA files were searched with DIA-NN 1.8.1 in library-free mode against
TriTrypDB release 51 TREU927 proteins plus the Tb427 expression-site
proteins and the GFP-BLA marker (parameters in `upstream/diann_search.md`).
Protein-group quantities were computed with MaxLFQ (diann R package) from
precursors with q-value ≤ 0.01 and protein-group q-value ≤ 0.01, and the
number of distinct peptides per protein group was counted
(`upstream/proteomics_maxlfq.R`). These two tables are
`data/proteomics/protein_groups.txt.gz` and `peptides_counts.txt`.

### Differential abundance (`scripts/03`, `scripts/04`)

1. Protein groups with at least 2 peptides were kept.
2. The 9 cell-extract samples were taken, and proteins with no value in
   any of them were removed (5,972 proteins remain).
3. Loading normalisation: each sample was scaled so that its median equals
   the mean of the sample medians.
4. Missing values were imputed separately for each temperature. Proteins
   missing in all 3 replicates of a temperature were given random values
   from a normal distribution truncated at 0, centred on the detection
   limit (median of the 10 smallest values of each sample), with standard
   deviation 20% of the limit. Other missing values were estimated with
   scikit-learn `IterativeImputer` and raised to the detection limit if
   below it.
5. Values were log2-transformed. limma fitted one mean per temperature, and
   the contrasts were tested with empirical Bayes moderation.
   P values were adjusted with the Benjamini-Hochberg method.

The published table is `data/proteomics/all_cellular_comparisons.csv.gz`.
limma reports `C34_C40_logFC` as 34 °C minus 40 °C; the figure script
reverses the sign so that positive means higher at 40 °C.

## RNA-seq

### Alignment and counting (`upstream/`)

Paired-end reads were trimmed with fastp 0.20.0 and aligned with bowtie2
(`--very-sensitive-local`) to two references: the TriTrypDB release 68
TREU927 genome, extended with the kinetoplast maxicircle, Tb427 VSG genes
and selectable marker genes (all alignments kept); and the Tb427
bloodstream-form expression sites (MAPQ ≥ 2). Duplicates were marked with
Picard 2.22.3 and kept. Read pairs were counted on CDS features per gene
with featureCounts 1.6.4 (`-p -B -C -M -O`). Commands:
`upstream/rnaseq_align_count.sh`.

### Count matrix (`scripts/01`)

The two count tables were stacked. The GFP-BLA marker and the 14
full-length expression-site VSG genes were removed, because their reads
are counted on the truncated VSG models. The counts of *AGO1*
(Tb927.10.10850) were set to the observed values in wild-type samples and
to 1 in all *ago1* null samples, where the gene is deleted. Genes without a
GC content and length from `bedtools nuc` were dropped
(10,131 genes remain).

### Differential expression (`scripts/02`)

1. edgeR `filterByExpr` (minimum count 10, minimum total count 30) kept
   expressed genes (9,396).
2. cqn computed offsets that correct for GC content and gene length.
3. A quasi-likelihood negative binomial GLM was fitted on all 27 samples,
   with 6 groups: WT and *ago1* (both lines pooled) at each temperature.
   Common, trended and tagwise dispersions were estimated.
4. The contrast WT 40 °C − WT 34 °C was tested (`glmQLFTest`), and P values
   were adjusted with the Benjamini-Hochberg method.

The published table is `data/rnaseq/WT40_vs_WT34.txt`.

## Gene sequences (`scripts/05`)

The annotation `data/annotation/merged_fixed_6.gtf.gz` combines the
TriTrypDB release 68 TREU927 gene models with transcript models curated by
hand in the Apollo genome browser, which extend genes into their UTRs. For
genes with a curated model, the TriTrypDB model was removed. From this
annotation and the TriTrypDB release 68 genome, the 5′UTR, CDS and 3′UTR
were extracted for genes with:

- exactly one CDS and one transcript;
- a transcript that extends more than 20 nt beyond the CDS at both ends;
- a CDS that starts with ATG.

6,812 genes pass (`data/sequences/`).

## Embeddings (`scripts/06`)

Each sequence was encoded with HyenaDNA medium-450k
(`LongSafari/hyenadna-medium-450k-seqlen-hf`, revision 42dedd4d374e).
Sequences were padded to 31,585 tokens, and the last hidden states were
averaged over the non-padding positions, giving 256 values per sequence
(`data/embeddings/`).

## Up and Down genes (`scripts/07`)

Both gene lists were restricted to the 6,923 genes of the
translation-efficiency data set of github.com/mtinti/decoding-gene-expression
(`data/gene_universe/`), and then to genes with embeddings.

| Data | Expression filter | Significance | Up | Down |
|---|---|---|---|---|
| Proteomics | average log2 abundance > 18 | adjusted P < 0.01 | logFC > 0.3 | logFC < −0.55 |
| RNA-seq | logCPM > 1 | FDR < 0.01 | logFC > 0.4 | logFC < −0.4 |

The proteomics thresholds are asymmetric so that the two groups have a
similar size. Final groups: proteomics 311 Up and 283 Down; RNA-seq 269 Up
and 242 Down.

## Prediction and ROC (`scripts/07`)

For each data set and each feature set (5′UTR + 3′UTR, 3′UTR only, 5′UTR
only), a random forest classifier (scikit-learn 1.6.1; 50 trees, maximum
depth 10, at least 2 samples per leaf, balanced class weights, random state
42) was evaluated by stratified 3-fold cross-validation with shuffling
(random state 42). The out-of-fold probability of Up was used to draw the
ROC curve and compute the AUC. To show the variability due to the random
split, the cross-validation was repeated with random states 0 to 19, and
the mean and standard deviation of the AUC are reported.

## Software

Docker image (`Dockerfile`): Python 3.12 with scikit-learn 1.6.1,
numpy 2.0.2, pandas 2.2.2, matplotlib 3.10.0, scipy 1.14.1 and biopython
1.85; R 4.4.3 with edgeR 4.4.0 (`glmQLFit(legacy = TRUE)`, the edgeR 3
method used for the published table), limma 3.62.1 and cqn 1.52.0.
Embedding image (`docker/Dockerfile.embeddings`): PyTorch 2.5.1 and
transformers 4.48.3.
