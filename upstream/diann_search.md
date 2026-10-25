# DIA-NN search

The proteomics quantities in `data/proteomics/` come from one DIA-NN 1.8.1
run on 18 DIA files (cell extract and culture supernatant at 34, 37 and
40 °C, 3 replicates each). The figure uses the 9 cell-extract files.

## Search database

`TriTrypDB-51_TbruceiTREU927_AnnotatedProteins_added_BES_GFP-BLA.fasta`:
TriTrypDB release 51 *T. brucei* TREU927 annotated proteins, plus the Tb427
bloodstream-form expression-site proteins and the GFP-BLA marker
(11,363 entries, MD5 `1007d6cb2905b8b0f6591490cd3f1df3`). It is deposited
with the raw data.

## Command

```
diann.exe --f <each of the 18 .raw.dia files> --lib --threads 32 --verbose 1
  --out report.tsv --qvalue 0.01 --matrices --out-lib report-lib.tsv
  --gen-spec-lib --predictor
  --fasta TriTrypDB-51_TbruceiTREU927_AnnotatedProteins_added_BES_GFP-BLA.fasta
  --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision
  --cut K*,R* --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30
  --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4
  --unimod4 --individual-mass-acc --individual-windows --reanalyse
  --relaxed-prot-inf --smart-profiling --peak-center --no-ifs-removal
```

In words: library-free search with an in-silico predicted library, trypsin
with at most 1 missed cleavage, peptides of 7 to 30 residues, precursor
charge 1 to 4, fixed carbamidomethylation of cysteine, N-terminal
methionine excision, 1% precursor FDR, and match-between-runs
(`--reanalyse`).

## From the report to protein quantities

`upstream/proteomics_maxlfq.R` turns `report.tsv` into
`protein_groups.txt` (MaxLFQ) and `peptides_counts.txt`.
