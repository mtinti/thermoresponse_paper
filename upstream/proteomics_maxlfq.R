# Protein quantities and peptide counts from the DIA-NN report.
#
# - protein_groups.txt: MaxLFQ per protein group from Precursor.Normalised,
#   precursor and protein-group q-value <= 0.01.
# - peptides_counts.txt: number of distinct peptides per protein group,
#   filtered on Global.Q.Value <= 0.01 and PG.Q.Value <= 0.01.
#
# Needs the diann R package (github.com/vdemichev/diann-rpackage) and the
# DIA-NN report.tsv from the search described in upstream/diann_search.md. This is the code that
# produced data/proteomics/; it is documented here and not part of the
# Docker run.
#
# Usage:
#   Rscript upstream/proteomics_maxlfq.R REPORT_TSV OUT_DIR

library(diann)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript upstream/proteomics_maxlfq.R REPORT_TSV OUT_DIR")
dir.create(args[2], showWarnings = FALSE, recursive = TRUE)

df <- diann_load(args[1])
passing <- df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01, ]

protein.groups <- diann_maxlfq(passing, group.header = "Protein.Group",
                               id.header = "Precursor.Id",
                               quantity.header = "Precursor.Normalised")
write.table(protein.groups, file = file.path(args[2], "protein_groups.txt"), sep = "\t")

pep <- unique(df[df$Global.Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,
                 c("Stripped.Sequence", "Protein.Group")])
write.table(table(pep$Protein.Group), file = file.path(args[2], "peptides_counts.txt"), sep = "\t")
