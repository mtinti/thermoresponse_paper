# Reproduce the analysis. Run inside the Docker image (see README.md):
#
#   make all          every step below, in order
#   make check-data   checksums of the input files
#   make rnaseq       RNA-seq count matrix and differential expression
#   make proteomics   proteomics imputation and differential abundance
#   make sequences    5'UTR, CDS and 3'UTR sequences from annotation and genome
#   make figure       ROC figure from the published intermediate tables
#   make figure-rebuilt   ROC figure from the tables rebuilt by this run
#   make embeddings   HyenaDNA embeddings (GPU image only, not part of all)

PYTHON ?= python
RSCRIPT ?= Rscript
OUT ?= results
REPEATS ?= 20

.PHONY: all check-data rnaseq proteomics sequences figure figure-rebuilt embeddings clean

all: check-data rnaseq proteomics sequences figure figure-rebuilt

check-data:
	$(PYTHON) scripts/verify.py inputs

rnaseq:
	mkdir -p $(OUT)/rnaseq
	$(PYTHON) scripts/01_rnaseq_matrix.py data/rnaseq $(OUT)/rnaseq/counts_with_gc.csv
	$(RSCRIPT) scripts/02_rnaseq_de.R $(OUT)/rnaseq/counts_with_gc.csv $(OUT)/rnaseq
	$(PYTHON) scripts/verify.py rnaseq $(OUT)/rnaseq/WT40_vs_WT34.txt

proteomics:
	$(PYTHON) scripts/03_proteomics_prepare.py data/proteomics $(OUT)/proteomics
	$(RSCRIPT) scripts/04_proteomics_de.R $(OUT)/proteomics $(OUT)/proteomics/all_cellular_comparisons.csv
	$(PYTHON) scripts/verify.py proteomics $(OUT)/proteomics/all_cellular_comparisons.csv

sequences:
	mkdir -p $(OUT)/sequences
	gzip -dc data/annotation/merged_fixed_6.gtf.gz > $(OUT)/sequences/merged_fixed_6.gtf
	gzip -dc data/reference/TriTrypDB-68_TbruceiTREU927_Genome.fasta.gz > $(OUT)/sequences/genome.fasta
	$(PYTHON) scripts/05_extract_sequences.py $(OUT)/sequences/merged_fixed_6.gtf $(OUT)/sequences/genome.fasta $(OUT)/sequences
	$(PYTHON) scripts/verify.py sequences $(OUT)/sequences

figure:
	$(PYTHON) scripts/07_figure_roc.py --out $(OUT)/figure --repeats $(REPEATS)
	$(PYTHON) scripts/verify.py figure $(OUT)/figure

figure-rebuilt: rnaseq proteomics
	$(PYTHON) scripts/07_figure_roc.py --out $(OUT)/figure_rebuilt --repeats $(REPEATS) \
		--proteomics $(OUT)/proteomics/all_cellular_comparisons.csv \
		--transcriptomics $(OUT)/rnaseq/WT40_vs_WT34.txt

embeddings:
	mkdir -p $(OUT)/embeddings
	for region in 5utr cds 3utr; do \
		$(PYTHON) scripts/06_compute_embeddings.py data/sequences/trypanosoma_sequences_$$region.csv.gz \
			$(OUT)/embeddings/embeddings_927_hyenadna-medium-450k-seqlen-hf_trypanosoma_sequences_$$region.csv.npy; \
	done

clean:
	rm -rf $(OUT)
