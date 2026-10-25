# Differential expression of the RNA-seq data with edgeR and cqn.
#
# filterByExpr removes lowly expressed genes; cqn corrects for gene GC
# content and length; a quasi-likelihood GLM is fitted on all 27 samples
# (WT and the two ago1 null lines at 34, 37 and 40 C) and the contrast
# WT40 - WT34 is tested. The two ago1 lines share the groups MUT34, MUT37
# and MUT40, and all samples inform the dispersion estimates.
#
# The published table was made with edgeR 3.x. edgeR 4 changed the default
# quasi-likelihood method, so the fit uses legacy = TRUE where available.
# cqn draws random numbers; the script sets a seed so that its output is
# stable. Without a seed, reruns differ in a few genes at the logFC cutoff.
#
# Usage:
#   Rscript scripts/02_rnaseq_de.R COUNTS_CSV OUT_DIR [SEED]

suppressPackageStartupMessages({
  library(edgeR)
  library(cqn)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2, 3)) stop("Usage: Rscript scripts/02_rnaseq_de.R COUNTS_CSV OUT_DIR [SEED]")
set.seed(if (length(args) == 3) as.integer(args[3]) else 1)
dir.create(args[2], showWarnings = FALSE, recursive = TRUE)

d <- read.csv(args[1], row.names = 1, check.names = FALSE)
gc <- d[, c("length", "gccontent")]
counts <- as.matrix(d[, setdiff(colnames(d), c("length", "gccontent"))])

group <- factor(c(rep(c("WT34", "WT37", "WT40"), each = 3),
                  rep(rep(c("MUT34", "MUT37", "MUT40"), each = 3), 2)))

y <- DGEList(counts = counts, group = group)
keep <- filterByExpr(y, min.count = 10, min.total.count = 30)
counts <- y$counts[keep, ]
gc <- gc[rownames(counts), ]
size_factors <- colSums(counts)

cq <- cqn(counts, lengths = gc$length, x = gc$gccontent,
          sizeFactors = size_factors, verbose = FALSE)

design <- model.matrix(~0 + group)
y <- DGEList(counts = counts, lib.size = size_factors, group = group)
y$offset <- cq$glm.offset
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- if ("legacy" %in% names(formals(glmQLFit.DGEList))) {
  glmQLFit(y, design, legacy = TRUE)
} else {
  glmQLFit(y, design)
}

contrast <- makeContrasts(groupWT40 - groupWT34, levels = design)
table <- topTags(glmQLFTest(fit, contrast = contrast),
                 n = Inf, sort.by = "none", adjust.method = "BH")$table
out <- file.path(args[2], "WT40_vs_WT34.txt")
write.table(table, file = out, sep = "\t", quote = FALSE, row.names = TRUE)
cat(out, ":", nrow(table), "genes. edgeR", as.character(packageVersion("edgeR")),
    "cqn", as.character(packageVersion("cqn")), "\n")
