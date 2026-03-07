############################################################
# Script: 01_remove_unwanted_variation.R
# Purpose: Remove unwanted technical variation from RNA-seq
#          count data using the RUVSeq framework.
# Input: salmon.merged.gene_counts.tsv
#        samplesheet.tsv
# Output: count matrix corrected
# Author: Bruna Toledo
############################################################


############################################################
# 1. Load libraries
############################################################

suppressPackageStartupMessages({
  library(RUVSeq)
  library(EDASeq)
  library(edgeR)
  library(RColorBrewer)
  library(readr)
  library(dplyr)
})


############################################################
# 2. Define input/output paths
############################################################

salmon_counts <- "salmon.merged.gene_counts.tsv"
samples_tsv   <- "samplesheet.tsv"

out_dir <- "1_RUVseq_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_base   <- "salmon.merged.gene_counts_RUVseq"
out_prefix <- file.path(out_dir, out_base)


############################################################
# 3. Load count matrix
############################################################

counts_raw <- read_tsv(salmon_counts, show_col_types = FALSE) %>%
  as.data.frame()

rownames(counts_raw) <- counts_raw$gene_id
counts_raw <- counts_raw[,3:14]
counts_raw <- as.matrix(counts_raw)


############################################################
# 4. Load and align metadata
############################################################

metadata <- read_tsv(samples_tsv, show_col_types = FALSE) %>%
  mutate(
    subject   = factor(subject),
    condition = factor(condition),
    batch     = as.character(batch)
  )

sample_ids <- as.character(metadata$sample)

stopifnot(setequal(colnames(counts_raw), sample_ids))

counts_raw <- counts_raw[, sample_ids]

stopifnot(identical(colnames(counts_raw), sample_ids))


############################################################
# 5. Filter low-abundance genes
############################################################

keep_abundance <- rowSums(counts_raw >= 25, na.rm = TRUE) >= 3
keep_nonNA     <- rowSums(!is.na(counts_raw)) >= (ncol(counts_raw) * 0.5)

keep_both <- keep_abundance & keep_nonNA

counts_ftd <- counts_raw[keep_both, , drop = FALSE]
counts_ftd <- round(counts_ftd)


############################################################
# 6. Create SeqExpressionSet
############################################################

set_raw <- newSeqExpressionSet(
  counts = counts_ftd,
  phenoData = data.frame(
    row.names = metadata$sample,
    condition = metadata$condition,
    subject   = metadata$subject,
    batch     = metadata$batch
  )
)


############################################################
# 7. Exploratory PCA (raw counts)
############################################################

colors <- brewer.pal(n = nlevels(metadata$condition), name = "Set2")

png(paste0(out_prefix, "-PCA_raw.png"),
    width = 1200, height = 1000, res = 150)

plotPCA(set_raw, col = colors[metadata$condition], cex = 1.2)

dev.off()


############################################################
# 8. Between-lane normalization (upper quartile)
############################################################

set_uq <- betweenLaneNormalization(set_raw, which = "upper")

png(paste0(out_prefix, "-PCA_uq.png"),
    width = 1200, height = 1000, res = 150)

plotPCA(set_uq, col = colors[metadata$condition], cex = 1.2)

dev.off()


############################################################
# 9. Remove unwanted variation (RUVr)
############################################################

design <- model.matrix(~ condition, data = pData(set_uq))

y <- DGEList(
  counts = EDASeq::counts(set_uq),
  group  = metadata$condition
)

y <- calcNormFactors(y, method = "upperquartile")

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)

resids <- residuals(fit, type = "deviance")


############################################################
# 10. Estimate RUV factors
############################################################

k <- 2

controls <- rownames(EDASeq::counts(set_uq))

set_ruvr <- RUVr(
  set_uq,
  controls,
  k = k,
  res = resids
)

pData(set_ruvr)


############################################################
# 11. PCA after RUV correction
############################################################

png(paste0(out_prefix, "-PCA_ruvr_k2.png"),
    width = 1200, height = 1000, res = 150)

plotPCA(set_ruvr, col = colors[metadata$condition], cex = 1.2)

dev.off()


############################################################
# 12. Save workspace
############################################################

save.image(file = "RUVseq_workspace_f253-k2.RData")


############################################################
# End of script
############################################################