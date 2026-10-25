# Differential protein abundance between temperatures with limma.
#
# Fits one mean per temperature (34, 37, 40 C) to the log2 imputed matrix
# from 03_proteomics_prepare.py and tests the contrasts C34-C37, C34-C40 and
# C37-C40 with empirical Bayes moderation. The output has the same layout as
# data/proteomics/all_cellular_comparisons.csv.gz: the protein quantities
# of all 18 samples, then logFC, AveExpr, t, P.Value, adj.P.Val (Benjamini-
# Hochberg) and B for each contrast, prefixed with the contrast name.
# logFC is the first temperature minus the second, so C34_C40_logFC < 0
# means higher at 40 C.
#
# Usage:
#   Rscript scripts/04_proteomics_de.R PREPARED_DIR OUT_CSV

suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript scripts/04_proteomics_de.R PREPARED_DIR OUT_CSV")

x <- as.matrix(read.csv(file.path(args[1], "proteomics_log2_imputed.csv"),
                        row.names = 1, check.names = FALSE))
quant <- read.csv(file.path(args[1], "proteomics_quantities.csv"),
                  row.names = 1, check.names = FALSE)

temperature <- factor(sub("\\..*", "", colnames(x)), levels = c("C34", "C37", "C40"))
design <- model.matrix(~0 + temperature)
colnames(design) <- levels(temperature)
contrasts <- makeContrasts(C34_C37 = C34 - C37, C34_C40 = C34 - C40,
                           C37_C40 = C37 - C40, levels = design)
fit <- eBayes(contrasts.fit(lmFit(x, design), contrasts))

out <- quant
for (name in colnames(contrasts)) {
  table <- topTable(fit, coef = name, sort.by = "none", number = Inf)
  colnames(table) <- paste(name, colnames(table), sep = "_")
  out <- merge(out, table, by = "row.names", all = TRUE)
  rownames(out) <- out$Row.names
  out$Row.names <- NULL
}
write.csv(out, args[2])
cat(args[2], ":", sum(!is.na(out$C34_C40_logFC)), "proteins tested\n")
