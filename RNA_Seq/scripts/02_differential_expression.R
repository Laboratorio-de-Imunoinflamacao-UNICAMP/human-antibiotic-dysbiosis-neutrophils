############################################################
# Script: 02_differential_expression.R
# Purpose: Perform differential gene expression analysis
#          using DESeq2 on RUV-corrected RNA-seq counts.
#
# Input:  RUVseq_workspace_f253-k2.RData
# Output: results/tables/differential_expression/
# Author: Bruna Toledo
############################################################

############################################################
# 1. Load libraries
############################################################

suppressPackageStartupMessages({
  library(ashr)
  library(readr)
  library(dplyr)
  library(DESeq2)
  library(writexl)
  library(stringr)
  library(rtracklayer)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})


############################################################
# 2. Define input/output paths and parameters
############################################################

out_dir <- "1_DEG_tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ruv_workspace <- "RUVseq_workspace_f253-k2.RData"
tx2gene_file  <- "tx2gene.tsv"
gtf_file      <- "gencode.v49.primary_assembly.annotation.gtf"

id_col   <- "gene_id"
lfc_col  <- "log2FoldChange"
padj_col <- "padj"

alpha   <- 0.05
lfc_thr <- log2(1.5)


############################################################
# 3. Load input files
############################################################

load(ruv_workspace)

gene_anno <- read_tsv(tx2gene_file, show_col_types = FALSE) %>%
  as.data.frame() %>%
  mutate(
    transcript_id = as.character(transcript_id),
    gene_id       = as.character(gene_id),
    gene_name     = as.character(gene_name)
  )


############################################################
# 4. Define helper functions
############################################################

get_results <- function(dds, contrast, alpha = 0.05, shrink = TRUE) {
  res <- DESeq2::results(dds, contrast = contrast, alpha = alpha)
  
  if (shrink) {
    if (!requireNamespace("ashr", quietly = TRUE)) {
      stop("Package 'ashr' is required for LFC shrinkage.")
    }
    res <- DESeq2::lfcShrink(dds, contrast = contrast, type = "ashr")
  }
  
  out <- as.data.frame(res)
  out$gene_id <- rownames(out)
  out <- out[, c("gene_id", setdiff(names(out), "gene_id"))]
  
  out
}

filter_deg <- function(df,
                       lfc_col = "log2FoldChange",
                       padj_col = "padj",
                       alpha = 0.05,
                       lfc_thr = log2(1.5)) {
  stopifnot(all(c(lfc_col, padj_col) %in% colnames(df)))
  
  df %>%
    mutate(
      regulation = ifelse(.data[[lfc_col]] > 0, "Up", "Down")
    ) %>%
    filter(!is.na(.data[[padj_col]])) %>%
    filter(.data[[padj_col]] < alpha, abs(.data[[lfc_col]]) >= lfc_thr)
}

write_tsv2 <- function(df, path) {
  readr::write_tsv(df, path)
}


############################################################
# 5. Prepare DESeq2 colData
############################################################

cd <- as.data.frame(pData(set_ruvr), stringsAsFactors = FALSE)

cd$subject   <- factor(cd$subject)
cd$condition <- factor(cd$condition)
cd$condition <- relevel(cd$condition, ref = "Pre-Abx")

use_W <- c("W_2")

rhs <- c("subject", use_W, "condition")
design_formula <- as.formula(paste("~", paste(rhs, collapse = " + ")))


############################################################
# 6. Build DESeq2 object and run model
############################################################

cts <- round(EDASeq::counts(set_ruvr))

stopifnot(identical(colnames(cts), rownames(cd)))

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData   = cd,
  design    = design_formula
)

dds <- DESeq2::DESeq(
  dds,
  test    = "Wald",
  fitType = "parametric"
)


############################################################
# 7. Define contrasts
############################################################

contrasts <- list(
  PostAbx_vs_PreAbx       = c("condition", "Post-Abx",     "Pre-Abx"),
  PreAbxLPS_vs_PreAbx     = c("condition", "Pre-Abx+LPS",  "Pre-Abx"),
  PostAbxLPS_vs_PostAbx   = c("condition", "Post-Abx+LPS", "Post-Abx"),
  PostAbxLPS_vs_PreAbxLPS = c("condition", "Post-Abx+LPS", "Pre-Abx+LPS"),
  PostAbx_vs_PreAbxLPS    = c("condition", "Post-Abx",     "Pre-Abx+LPS")
)


############################################################
# 8. Compute differential expression results
############################################################

DEG_res_list <- lapply(
  contrasts,
  function(ct) get_results(dds, ct, alpha = alpha, shrink = TRUE)
)


############################################################
# 9. Load and prepare gene annotation
############################################################

gtf <- import(gtf_file)

gene_annot <- gtf %>%
  as.data.frame() %>%
  filter(type == "gene") %>%
  transmute(
    gene_id   = gene_id,
    gene_name = gene_name,
    gene_type = gene_type
  ) %>%
  distinct()

gene_map <- setNames(
  as.character(gene_annot$gene_name),
  as.character(gene_annot$gene_id)
)

gene_type_map <- setNames(
  gene_annot$gene_type,
  gene_annot$gene_id
)


############################################################
# 10. Annotate DEG tables
############################################################

DEG_res_list <- lapply(DEG_res_list, function(df) {
  stopifnot("gene_id" %in% names(df))
  
  df$gene_name <- gene_map[as.character(df$gene_id)]
  df$gene_type <- gene_type_map[as.character(df$gene_id)]
  
  df
})


############################################################
# 11. Export full DEG tables
############################################################

invisible(
  mapply(function(name, df) {
    write_tsv2(df, file.path(out_dir, paste0(name, ".deseq2.results.tsv")))
  }, names(DEG_res_list), DEG_res_list)
)


############################################################
# 12. Filter significant DEGs
############################################################

deg_list_sig <- lapply(DEG_res_list, function(df) {
  filter_deg(
    df,
    lfc_col  = lfc_col,
    padj_col = padj_col,
    alpha    = alpha,
    lfc_thr  = lfc_thr
  )
})


############################################################
# 13. Export filtered DEG tables
############################################################

invisible(
  mapply(function(name, df) {
    write_tsv2(df, file.path(out_dir, paste0(name, ".deseq2.results_sig.tsv")))
  }, names(deg_list_sig), deg_list_sig)
)

write_xlsx(
  DEG_res_list,
  path = file.path(out_dir, "DEG_all.xlsx")
)

write_xlsx(
  deg_list_sig,
  path = file.path(out_dir, "DEG_filtered.xlsx")
)


############################################################
# 14. Export normalized counts
############################################################

norm_counts <- as.data.frame(set_ruvr@assayData[["normalizedCounts"]])

norm_counts_anno <- norm_counts
norm_counts_anno$gene_id   <- rownames(norm_counts)
norm_counts_anno$gene_name <- gene_map[norm_counts_anno$gene_id]

norm_counts_anno <- norm_counts_anno[, c(
  "gene_id",
  "gene_name",
  setdiff(names(norm_counts_anno), c("gene_id", "gene_name"))
)]

readr::write_tsv(
  norm_counts_anno,
  file.path(out_dir, "RUV-f253-k2_normalised-counts_annotated.tsv")
)


############################################################
# 15. Save objects for downstream analyses
############################################################

save(
  DEG_res_list,
  deg_list_sig,
  file = file.path(out_dir, "DEG-lists_f253-k2.RData")
)

save(
  norm_counts,
  norm_counts_anno,
  counts_raw,
  counts_ftd,
  file = file.path(out_dir, "Counts_f253-k2.RData")
)

save(
  metadata,
  gene_anno,
  file = file.path(out_dir, "Metadata.RData")
)


############################################################
# 16. Save full workspace
############################################################

save.image(file = file.path(out_dir, "DEG-workspace_f253-k2.RData"))



############################################################
# End of script
############################################################