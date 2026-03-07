############################################################
# Script: 03_deg_visualization_and_enrichment.R
# Purpose: Generate RNA-seq QC plots, DEG visualizations,
#          and functional enrichment analyses.
# Input: 1_DEG_tables/
# Output: QC plots, DEG summaries, enrichment results
# Author: Bruna Toledo
############################################################

############################################################
# 1. Load libraries
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(forcats)
  library(writexl)
  library(stringr)
  library(openxlsx)
  library(pheatmap)
  library(ReactomePA)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(RColorBrewer)
  library(clusterProfiler)
  library(EnhancedVolcano)
  library(AnnotationDbi)
  library(readr)
  library(tibble)
  library(matrixStats)
  library(ggplotify)
})

############################################################
# 2. Define input/output paths
############################################################

deg_tables_dir   <- "1_DEG_tables"
qc_dir           <- "2_QC_output"
deg_overview_dir <- "3_DEG_overview"
enrich_dir       <- "4_DEG_Enrichment"
tf_dir           <- "5_TFs"

counts_file   <- file.path(deg_tables_dir, "Counts_f253-k2.RData")
deg_lists_file <- file.path(deg_tables_dir, "DEG-lists_f253-k2.RData")
metadata_file <- file.path(deg_tables_dir, "Metadata.RData")

dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(deg_overview_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cluster_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tf_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
# 3. Load input data
############################################################

load(counts_file)
load(deg_lists_file)
load(metadata_file)

############################################################
# 4. Helper functions
############################################################

stop_if_non_numeric <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) storage.mode(x) <- "numeric"
  x
}

ggsave2 <- function(p, file, width = 6.5, height = 5, out_dir = ".") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  stem <- tools::file_path_sans_ext(basename(file))
  path_tiff <- file.path(out_dir, paste0(stem, ".tiff"))
  path_svg  <- file.path(out_dir, paste0(stem, ".svg"))
  
  ggplot2::ggsave(
    filename = path_tiff, plot = p,
    width = width, height = height, dpi = 300,
    device = "tiff", compression = "lzw"
  )
  ggplot2::ggsave(
    filename = path_svg, plot = p,
    width = width, height = height, dpi = 300,
    device = "svg"
  )
  invisible(list(tiff = path_tiff, svg = path_svg))
}

log2p1 <- function(m) log2(m + 1)

map_aes <- function(meta) {
  aes_list <- list(color = NULL, shape = NULL)
  if ("condition" %in% names(meta)) aes_list$color <- meta$condition
  if ("batch" %in% names(meta))     aes_list$shape <- meta$batch
  aes_list
}

do_pca_plot <- function(m, meta, title) {
  m <- stop_if_non_numeric(m)
  X <- t(m)
  keep <- apply(X, 2, sd, na.rm = TRUE) > 0
  X <- X[, keep, drop = FALSE]
  pc <- prcomp(X, center = TRUE, scale. = FALSE)
  pc_df <- as.data.frame(pc$x[, 1:2, drop = FALSE])
  pc_df$sample <- rownames(pc$x)
  df <- left_join(pc_df, meta, by = c("sample" = "sample"))
  
  var_exp <- (pc$sdev^2) / sum(pc$sdev^2)
  xlab <- sprintf("PC1 (%.1f%%)", 100 * var_exp[1])
  ylab <- sprintf("PC2 (%.1f%%)", 100 * var_exp[2])
  
  aes_map <- map_aes(df)
  p <- ggplot(df, aes(PC1, PC2)) +
    {
      if (!is.null(aes_map$color) && !is.null(aes_map$shape)) {
        geom_point(aes(color = aes_map$color, shape = aes_map$shape), size = 3)
      } else if (!is.null(aes_map$color)) {
        geom_point(aes(color = aes_map$color), size = 3)
      } else {
        geom_point(size = 3)
      }
    } +
    labs(title = title, x = xlab, y = ylab, color = "condition", shape = "batch") +
    theme_bw()
  p
}

do_mds_plot <- function(m, meta, title) {
  m <- stop_if_non_numeric(m)
  X <- t(m)
  keep <- apply(X, 2, sd, na.rm = TRUE) > 0
  X <- X[, keep, drop = FALSE]
  d <- dist(scale(X))
  mds <- cmdscale(d, k = 2)
  df <- data.frame(MDS1 = mds[, 1], MDS2 = mds[, 2], sample = rownames(mds))
  df <- left_join(df, meta, by = "sample")
  
  aes_map <- map_aes(df)
  p <- ggplot(df, aes(MDS1, MDS2)) +
    {
      if (!is.null(aes_map$color) && !is.null(aes_map$shape)) {
        geom_point(aes(color = aes_map$color, shape = aes_map$shape), size = 3)
      } else if (!is.null(aes_map$color)) {
        geom_point(aes(color = aes_map$color), size = 3)
      } else {
        geom_point(size = 3)
      }
    } +
    labs(title = title, color = "condition", shape = "batch") +
    theme_bw()
  p
}

long_expr <- function(m, meta) {
  m <- stop_if_non_numeric(m)
  df <- as.data.frame(m, check.names = FALSE)
  df$gene <- rownames(df)
  df <- pivot_longer(df, -gene, names_to = "sample", values_to = "expr")
  left_join(df, meta, by = "sample")
}

.check_cols <- function(df, cols) stopifnot(all(cols %in% names(df)))

alpha   <- 0.05
lfc_thr <- log2(1.5)

ensure_regulation <- function(df) {
  if (!"regulation" %in% names(df)) {
    stopifnot(all(c("log2FoldChange", "padj") %in% names(df)))
    df <- df %>%
      mutate(
        regulation = case_when(
          padj < alpha & log2FoldChange >=  lfc_thr ~ "Up",
          padj < alpha & log2FoldChange <= -lfc_thr ~ "Down",
          TRUE ~ NA_character_
        )
      )
  }
  df %>% filter(!is.na(regulation))
}

count_up_down <- function(df, dataset_name) {
  df <- ensure_regulation(df)
  tibble(
    Dataset = dataset_name,
    Status  = df$regulation
  ) |>
    count(Dataset, Status, name = "Count")
}

nm_pretty <- function(x) {
  x <- as.character(x)
  if (exists("label_map")) {
    lab <- unname(label_map[x])
    lab[is.na(lab) | lab == ""] <- x[is.na(lab) | lab == ""]
    return(lab)
  }
  x
}

ev_volcano <- function(df,
                       title,
                       select_labels = character(0),
                       label_col = "gene_id",
                       p_col = "padj",
                       fc_col = "log2FoldChange",
                       p_cut = alpha,
                       fc_cut = lfc_thr,
                       xlim_override = NULL) {
  stopifnot(all(c(label_col, p_col, fc_col) %in% names(df)))
  
  lab_plain <- as.character(df[[label_col]])
  pvals     <- as.numeric(df[[p_col]])
  lfcs      <- as.numeric(df[[fc_col]])
  
  ok   <- !is.na(pvals) & !is.na(lfcs)
  up_n <- sum(ok & pvals < p_cut & lfcs >=  fc_cut)
  dn_n <- sum(ok & pvals < p_cut & lfcs <= -fc_cut)
  cap_txt <- sprintf("DEG: Up = %d, Down = %d", up_n, dn_n)
  
  keyvals <- ifelse(lfcs < -fc_cut & pvals < p_cut, "royalblue",
                    ifelse(lfcs >  fc_cut & pvals < p_cut, "red2", "grey30"))
  keyvals[is.na(keyvals)] <- "grey30"
  names(keyvals)[keyvals == "red2"]      <- "Upregulated"
  names(keyvals)[keyvals == "royalblue"] <- "Downregulated"
  names(keyvals)[keyvals == "grey30"]    <- "Not significant"
  
  sel    <- unique(select_labels)
  is_sel <- lab_plain %in% sel
  
  lab_disp <- lab_plain
  if (any(is_sel)) lab_disp[is_sel] <- paste0("italic('", lab_plain[is_sel], "')")
  
  df2 <- df
  df2[[label_col]] <- lab_plain
  df2[[p_col]]     <- pvals
  df2[[fc_col]]    <- lfcs
  
  xlim_final <- if (!is.null(xlim_override)) xlim_override else c(-5, 5)
  
  EnhancedVolcano::EnhancedVolcano(
    toptable = df2,
    lab = lab_disp,
    selectLab = lab_disp[is_sel],
    parseLabels = TRUE,
    x = fc_col, y = p_col,
    pCutoff = p_cut, FCcutoff = fc_cut,
    xlim = xlim_final,
    colCustom = keyvals,
    drawConnectors = TRUE, widthConnectors = 0.5, colConnectors = "grey40",
    max.overlaps = Inf, boxedLabels = TRUE, colAlpha = 0.3,
    title = paste(title),
    subtitle = "Differential gene expression - FC cutoff: 1.5; p-value cutoff: 0.05",
    caption = cap_txt,
    xlab = bquote(Log[2] ~ "fold change"),
    ylab = bquote(-Log[10] ~ adjusted ~ italic(P)),
    titleLabSize = 8, axisLabSize = 8, legendLabSize = 7,
    subtitleLabSize = 7, captionLabSize = 7, pointSize = 2, labSize = 3,
    gridlines.major = FALSE, gridlines.minor = FALSE,
    legendPosition = "top", border = "partial"
  )
}

compact_label <- function(x) {
  x <- gsub("Post-Abx \\+ LPS", "PostLPS", x, fixed = TRUE)
  x <- gsub("Pre-Abx \\+ LPS",  "PreLPS",  x, fixed = TRUE)
  x <- gsub("Post-Abx",         "Post",    x, fixed = TRUE)
  x <- gsub("Pre-Abx",          "Pre",     x, fixed = TRUE)
  x <- gsub(" vs ",             "_vs_",    x, fixed = TRUE)
  x <- gsub("[^A-Za-z0-9_]+",   "",        x)
  substr(x, 1, 31)
}

sheet_name_pair <- function(a, b) {
  nm <- paste0(compact_label(a), "_AND_", compact_label(b))
  substr(nm, 1, 40)
}

pairwise_sheets <- function(sets_named) {
  stopifnot(length(sets_named) >= 2, !is.null(names(sets_named)))
  nms <- names(sets_named)
  cmb <- utils::combn(nms, 2, simplify = FALSE)
  sheets <- lapply(cmb, function(nm_pair) {
    g <- intersect(sets_named[[nm_pair[1]]], sets_named[[nm_pair[2]]])
    data.frame(gene_id = sort(unique(as.character(g))), stringsAsFactors = FALSE)
  })
  names(sheets) <- vapply(cmb, function(nm_pair) {
    sheet_name_pair(nm_pair[1], nm_pair[2])
  }, character(1))
  sheets
}

exclusive_per_contrast <- function(sets_named) {
  stopifnot(length(sets_named) >= 2, !is.null(names(sets_named)))
  lapply(names(sets_named), function(nm) {
    others <- sets_named[names(sets_named) != nm]
    unique(setdiff(sets_named[[nm]], unique(unlist(others))))
  }) |>
    `names<-`(names(sets_named))
}

write_exclusive_xlsx <- function(sets_named, path) {
  ex_sets <- exclusive_per_contrast(sets_named)
  sheets <- lapply(names(ex_sets), function(nm) {
    data.frame(gene_id = sort(as.character(ex_sets[[nm]])), stringsAsFactors = FALSE)
  })
  names(sheets) <- vapply(names(ex_sets), compact_label, character(1))
  writexl::write_xlsx(sheets, path = path)
}

sym_to_entrez <- function(symbols) {
  sy <- unique(as.character(symbols))
  mp <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = sy,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  out <- unname(mp[sy])
  out <- as.character(out)
  unique(out[!is.na(out) & nzchar(out)])
}

sanitize_name <- function(x, maxlen = 100) {
  x <- str_replace_all(x, "[^A-Za-z0-9_\\-\\+ ]", "_")
  x <- str_squish(x)
  str_trunc(x, maxlen)
}

run_ora <- function(entrez_ids) {
  if (length(entrez_ids) == 0) return(list(GO = NULL, KEGG = NULL, REACT = NULL))
  list(
    GO = tryCatch(
      enrichGO(
        gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
        ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE
      ),
      error = function(e) NULL
    ),
    KEGG = tryCatch(
      enrichKEGG(
        gene = entrez_ids, organism = "hsa",
        pvalueCutoff = 0.05, pAdjustMethod = "BH"
      ),
      error = function(e) NULL
    ),
    REACT = tryCatch(
      enrichPathway(
        gene = entrez_ids, organism = "human",
        pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE
      ),
      error = function(e) NULL
    )
  )
}

write_or_placeholder <- function(eobj, xlsx_path, label) {
  if (!is.null(eobj) && nrow(as.data.frame(eobj)) > 0) {
    writexl::write_xlsx(as.data.frame(eobj), xlsx_path)
  } else {
    msg <- paste("No enrichment terms found for", label)
    writexl::write_xlsx(list(Results = data.frame(Message = msg)), xlsx_path)
  }
}

.parse_ratio <- function(df) {
  if (!"GeneRatio" %in% names(df)) return(rep(NA_real_, nrow(df)))
  num <- suppressWarnings(as.numeric(sub("/.*$", "", df$GeneRatio)))
  den <- suppressWarnings(as.numeric(sub("^.*/", "", df$GeneRatio)))
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}

go_dotplot_df <- function(df,
                          include_desc,
                          regex = TRUE,
                          top_n = Inf,
                          title = NULL,
                          out_file = NULL) {
  req <- c("Description", "p.adjust", "Count")
  stopifnot(all(req %in% names(df)))
  
  df <- df %>%
    mutate(
      GeneRatio_num = .parse_ratio(cur_data()),
      GeneRatio_num = ifelse(is.na(GeneRatio_num),
                             Count / max(Count, na.rm = TRUE),
                             GeneRatio_num)
    )
  
  if (regex) {
    pat <- paste0("(", paste(include_desc, collapse = ")|("), ")")
    dff <- df %>% filter(str_detect(Description, regex(pat, ignore_case = TRUE)))
  } else {
    dff <- df %>% filter(Description %in% include_desc)
  }
  
  if (!is.finite(top_n)) top_n <- nrow(dff)
  dff <- dff %>% slice_head(n = top_n)
  if (nrow(dff) == 0) stop("No GO terms matched.")
  
  dff <- dff %>% arrange(p.adjust, desc(GeneRatio_num))
  dff$Description <- factor(dff$Description, levels = rev(unique(dff$Description)))
  
  p <- ggplot(dff, aes(x = GeneRatio_num, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust, alpha = 0.4)) +
    scale_size(range = c(0.5, 3)) +
    scale_x_continuous(name = "Gene ratio") +
    scale_y_discrete(name = NULL) +
    scale_color_gradient(name = "padj", low = "red2", high = "royalblue", trans = "reverse") +
    guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
    labs(title = title %||% "GO enrichment (selected terms)") +
    theme_bw(base_size = 6) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(family = "Arial"),
      plot.title = element_text(size = 6, family = "Arial"),
      axis.title = element_text(size = 6, family = "Arial"),
      axis.text = element_text(size = 6, family = "Arial"),
      legend.text = element_text(size = 6, family = "Arial")
    )
  
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    ggsave(paste0(out_file, ".tiff"), p, width = 8, height = 3, units = "cm",
           dpi = 300, device = "tiff", compression = "lzw")
    ggsave(paste0(out_file, ".svg"), p, width = 8, height = 3, units = "cm",
           dpi = 300, device = "svg")
  }
  
  p
}

run_cc <- function(symbol_sets, tag, out_dir) {
  cat("Running compareCluster for", tag, "...\n")
  
  entrez_sets <- lapply(symbol_sets, sym_to_entrez)
  entrez_sets <- entrez_sets[vapply(entrez_sets, length, integer(1)) > 0]
  
  if (length(entrez_sets) == 0) {
    cat("  - No valid gene sets for", tag, "(skipping)\n")
    return(invisible(NULL))
  }
  
  cc <- tryCatch(
    compareCluster(
      geneCluster = entrez_sets,
      fun = "enrichGO",
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    ),
    error = function(e) {
      message("  - Error in compareCluster for ", tag, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(cc)) return(invisible(NULL))
  
  cc_df <- as.data.frame(cc)
  
  if (nrow(cc_df) > 0) {
    p <- dotplot(cc, showCategory = 5) +
      scale_size(range = c(1, 5)) +
      ggtitle(paste("compareCluster GO:BP -", tag)) +
      theme_bw(base_size = 6) +
      theme(
        text = element_text(family = "Arial"),
        plot.title = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6)
      )
    
    base <- file.path(out_dir, paste0("compareCluster_GO_", tag))
    ggsave(paste0(base, ".svg"), p, width = 9, height = 8.5, units = "cm", dpi = 300, device = "svg")
    ggsave(paste0(base, ".tiff"), p, width = 9, height = 8.5, units = "cm", dpi = 300, device = "tiff")
    
    writexl::write_xlsx(cc_df, paste0(base, ".xlsx"))
    saveRDS(cc, paste0(base, ".rds"))
  } else {
    placeholder <- data.frame(Message = paste("No enriched terms found for", tag))
    writexl::write_xlsx(
      list(Results = placeholder),
      file.path(out_dir, paste0("compareCluster_GO_", tag, "_empty.xlsx"))
    )
    saveRDS(cc, file.path(out_dir, paste0("compareCluster_GO_", tag, "_empty.rds")))
  }
  
  invisible(cc)
}

save_plot <- function(plot, filename_base) {
  dir.create(dirname(filename_base), showWarnings = FALSE, recursive = TRUE)
  ggsave(paste0(filename_base, ".tiff"), plot, width = 10, height = 8, dpi = 300,
         device = "tiff", compression = "lzw")
  ggsave(paste0(filename_base, ".svg"), plot, width = 10, height = 8, dpi = 300,
         device = "svg")
}

perform_network_analysis <- function(gene_data, contrast_name, subset_label = "All",
                                     output_base = file.path(enrich_dir, "Plots", "Network_Results")) {
  clean_contrast <- sanitize_name(contrast_name)
  clean_subset   <- sanitize_name(subset_label)
  subfolder <- file.path(output_base, clean_contrast, clean_subset)
  dir.create(subfolder, showWarnings = FALSE, recursive = TRUE)
  
  sy <- if ("hgnc_symbol" %in% names(gene_data)) {
    gene_data$hgnc_symbol
  } else if ("gene_id" %in% names(gene_data)) {
    gene_data$gene_id
  } else {
    stop("No gene column found (expected 'gene_id' or 'hgnc_symbol').")
  }
  
  gene_list <- unique(as.character(sy[gene_data$padj < 0.05]))
  
  valid_sy <- tryCatch(keys(org.Hs.eg.db, keytype = "SYMBOL"), error = function(e) character(0))
  gene_list <- intersect(gene_list, valid_sy)
  
  if (!length(gene_list)) {
    message("No valid SYMBOLs for ", contrast_name, " - ", subset_label, " → Skipped.")
    return(NULL)
  }
  
  go_results <- tryCatch(
    enrichGO(
      gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
      ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
    net_plot <- cnetplot(go_results, showCategory = 10, circular = FALSE, colorEdge = TRUE) +
      ggtitle(paste("Gene-Concept Network -", contrast_name, "-", subset_label))
    save_plot(net_plot, file.path(subfolder, "Network"))
  }
  
  writexl::write_xlsx(
    if (!is.null(go_results)) as.data.frame(go_results) else data.frame(Message = "No enriched terms"),
    file.path(subfolder, paste0("GO_Network_", clean_contrast, "_", clean_subset, ".xlsx"))
  )
  
  go_results
}

scale_rows <- function(m) {
  m <- sweep(m, 1, rowMeans(m), FUN = "-")
  s <- matrixStats::rowSds(m)
  s[s == 0] <- 1
  sweep(m, 1, s, FUN = "/")
}

extract_deg_class <- function(datasets_filt, class_file, out_dir, class_label = "TF") {
  class_genes <- read.delim(class_file, stringsAsFactors = FALSE, header = TRUE) |>
    dplyr::pull(Symbol) |>
    as.character() |>
    unique() |>
    (\(x) x[!is.na(x) & nzchar(x)])()
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  class_results <- list()
  
  for (contrast in names(datasets_filt)) {
    message("Processing ", class_label, " DEGs for: ", contrast)
    
    deg <- datasets_filt[[contrast]] |>
      mutate(gene_id = as.character(gene_id)) |>
      filter(padj < 0.05, gene_id %in% class_genes)
    
    if (nrow(deg) == 0) {
      message("⚠️ No ", class_label, " DEGs found for ", contrast)
      next
    }
    
    deg_up   <- filter(deg, log2FoldChange > 0)
    deg_down <- filter(deg, log2FoldChange < 0)
    
    file_out <- file.path(out_dir, paste0(sanitize_name(contrast), "_", class_label, "_DEGs.xlsx"))
    write_xlsx(
      list(All = deg, Upregulated = deg_up, Downregulated = deg_down),
      path = file_out
    )
    
    class_results[[contrast]] <- list(
      All = deg,
      Upregulated = deg_up,
      Downregulated = deg_down
    )
  }
  
  invisible(class_results)
}

plot_deg_class <- function(df, contrast_name, class_label, out_dir = NULL,
                           color_up = "red2", color_down = "royalblue") {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  
  df <- df %>%
    mutate(
      regulation = case_when(
        log2FoldChange > 0 ~ "Up",
        log2FoldChange < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  for (dir in c("Up", "Down")) {
    sub_df <- df %>% filter(regulation == dir)
    if (nrow(sub_df) == 0) next
    
    sub_df <- sub_df %>%
      arrange(desc(abs(log2FoldChange))) %>%
      slice_head(n = 30)
    
    p <- ggplot(sub_df, aes(x = log10(baseMean), y = log2FoldChange)) +
      geom_point(color = ifelse(dir == "Up", color_up, color_down), size = 2, alpha = 0.3) +
      geom_text_repel(aes(label = gene_id), size = 3, max.overlaps = 15) +
      labs(
        title = paste(contrast_name, "-", class_label, "-", dir, "regulated"),
        x = "log10(Base Mean Expression)",
        y = "log2 Fold Change"
      ) +
      theme_minimal(base_size = 11)
    
    if (!is.null(out_dir)) {
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      file_base <- file.path(out_dir, paste0(sanitize_name(contrast_name), "_", class_label, "_", dir))
      ggsave(paste0(file_base, ".tiff"), plot = p, width = 6, height = 5, dpi = 300, compression = "lzw")
      ggsave(paste0(file_base, ".svg"), plot = p, width = 6, height = 5, dpi = 300)
    }
    
    print(p)
  }
}

sanitize_sheet_name <- function(x) {
  x <- str_replace_all(x, "[^a-zA-Z0-9_-]", "_")
  str_trunc(x, 31)
}

save_list_to_excel <- function(dataset_list, output_file = "DEG_Data.xlsx") {
  sheet_names <- sanitize_sheet_name(names(dataset_list))
  clean_list <- lapply(dataset_list, function(df) {
    if (is.data.frame(df)) df else data.frame(Message = "Not a data frame")
  })
  names(clean_list) <- sheet_names
  write_xlsx(clean_list, output_file)
  message("✅ Saved Excel file with ", length(clean_list), " sheets: ", output_file)
}

############################################################
# 5. Prepare numeric matrices
############################################################

counts_raw  <- stop_if_non_numeric(counts_raw)
counts_ftd  <- stop_if_non_numeric(counts_ftd)
norm_counts <- stop_if_non_numeric(norm_counts)

############################################################
# 6. Quality control plots
############################################################

## 6.1 Library sizes
libsizes <- data.frame(
  sample = colnames(counts_raw),
  raw    = colSums(counts_raw),
  ftd    = colSums(counts_ftd),
  norm   = colSums(norm_counts)
) %>%
  pivot_longer(-sample, names_to = "stage", values_to = "reads")

p_lib <- ggplot(libsizes, aes(x = sample, y = reads / 1e6, fill = stage)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Library sizes (millions)", x = "Sample", y = "Reads (M)", fill = "Stage") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave2(
  p_lib,
  "qc_library_sizes",
  width = 9,
  height = 4.5,
  out_dir = file.path(qc_dir, "1_Library-Sizes")
)

## 6.2 Boxplots
for (nm in c("Raw" = "counts_raw", "Filtered" = "counts_ftd", "Normalized" = "norm_counts")) {
  mat <- get(nm)
  df <- long_expr(log2p1(mat), metadata)
  
  p <- ggplot(df, aes(x = sample, y = expr)) +
    geom_boxplot(outlier.size = 0.2) +
    labs(
      title = paste("Boxplot log2(count+1) -", nm),
      x = "Sample", y = "log2(count+1)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave2(
    p,
    paste0("qc_boxplot_", nm),
    width = 9,
    height = 4.5,
    out_dir = file.path(qc_dir, "2_Boxplot")
  )
}

## 6.3 Density
df_den <- long_expr(log2p1(norm_counts), metadata)

p_den <- ggplot(df_den, aes(x = expr, group = sample)) +
  geom_density(alpha = 0.6) +
  labs(title = "Density of log2(count+1) - Normalized", x = "log2(count+1)", y = "density") +
  theme_bw()

ggsave2(
  p_den,
  "qc_density_normalized",
  width = 7,
  height = 5,
  out_dir = file.path(qc_dir, "3_Density-plot")
)

## 6.4 PCA and MDS
datasets <- list(
  Raw = counts_raw,
  Filtered = counts_ftd,
  Normalized = norm_counts
)
datasets_log <- lapply(datasets, log2p1)

plot_fns  <- list(PCA = do_pca_plot, MDS = do_mds_plot)
plot_dirs <- c(
  PCA = file.path(qc_dir, "5_PCA"),
  MDS = file.path(qc_dir, "4_MDS")
)

invisible(lapply(unique(plot_dirs), dir.create, showWarnings = FALSE, recursive = TRUE))

for (kind in names(plot_fns)) {
  fn <- plot_fns[[kind]]
  od <- plot_dirs[[kind]]
  
  for (nm in names(datasets_log)) {
    mat <- datasets_log[[nm]]
    title <- sprintf("%s - %s (log2+1)", kind, nm)
    p <- fn(mat, metadata, title)
    
    ggsave2(
      p,
      file = sprintf("qc_%s_%s", tolower(kind), tolower(nm)),
      out_dir = od
    )
  }
}

## 6.5 PCA normalized counts
pca <- prcomp(t(log2p1(norm_counts)), center = TRUE, scale. = FALSE)

ev <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 2)
pca_df <- transform(
  as.data.frame(pca$x[, 1:2]),
  sample    = rownames(pca$x),
  subject   = factor(metadata$subject),
  condition = factor(metadata$condition,
                     levels = c("Pre-Abx", "Pre-Abx+LPS", "Post-Abx", "Post-Abx+LPS"))
)

group_colors <- c(
  "Pre-Abx"      = "#804D80CC",
  "Pre-Abx+LPS"  = "#804D80CC",
  "Post-Abx"     = "#804D80CC",
  "Post-Abx+LPS" = "#804D80CC"
)

shape_map <- c(
  "Pre-Abx"      = 1,
  "Pre-Abx+LPS"  = 2,
  "Post-Abx"     = 19,
  "Post-Abx+LPS" = 17
)

pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = condition, label = subject)) +
  geom_point(size = 3) +
  geom_text(vjust = 2, size = 1) +
  scale_color_manual(values = group_colors, breaks = levels(pca_df$condition), name = NULL) +
  scale_shape_manual(values = shape_map, breaks = levels(pca_df$condition), name = NULL) +
  labs(
    title = "RNA-seq - PCA of Normalized Counts",
    x = sprintf("PC1 (%.2f%%)", ev[1]),
    y = sprintf("PC2 (%.2f%%)", ev[2])
  ) +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.2),
    axis.line.y.left   = element_line(color = "black", linewidth = 0.2),
    axis.ticks         = element_line(color = "black")
  )

ggsave(file.path(qc_dir, "PCA2.svg"), plot = pca_plot, width = 6.5, height = 6, units = "cm")

############################################################
# 7. Differential expression overview
############################################################

label_map <- c(
  PostAbx_vs_PreAbx           = "Post-Abx vs Pre-Abx",
  PreAbxLPS_vs_PreAbx         = "Pre-Abx + LPS vs Pre-Abx",
  PostAbxLPS_vs_PreAbxLPS     = "Post-Abx + LPS vs Pre-Abx + LPS",
  PostAbxLPS_vs_PostAbx       = "Post-Abx + LPS vs Post-Abx",
  PostAbx_vs_PreAbxLPS        = "Post-Abx vs Pre-Abx + LPS"
)

manual_labels <- list(
  PostAbx_vs_PreAbx = c("CCL2","CCL20","EDN1","CXCL8","IL1A","NFKBIA","NFKB1",
                        "IL1B","NR4A1","TNF","CSK","SLC16A6","FOLR3","CRYBG1","RHOQ"),
  PreAbxLPS_vs_PreAbx = c("PTGS2","EDN1","CXCL2","CCL2","NFKB1","IL1A","TNFAIP3",
                          "NR4A1","PTGER2","IL1B","NLRP3","CXCL3","CXCL8","TNF",
                          "SOD2","RIG1","LTBR","TRAF2"),
  PostAbxLPS_vs_PostAbx = c("F8A2"),
  PostAbxLPS_vs_PreAbxLPS = c("TNFRSF12A","SLA2","MS4A6A","IL9R","VNN3P","PRR5L","LAIR1"),
  PostAbx_vs_PreAbxLPS = c("F8A2","TNFRSF12A","VNN3P","IL9R","CNR2")
)

## 7.1 Volcano plots
volcano_dir <- file.path(deg_overview_dir, "Volcano_Plot")
dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)

invisible(lapply(names(DEG_res_list), function(nm) {
  df  <- DEG_res_list[[nm]]
  sel <- if (!is.null(manual_labels[[nm]])) manual_labels[[nm]] else character(0)
  
  xlim_use <- if (nm %in% c("PostAbxLPS_vs_PostAbx", "PostAbxLPS_vs_PreAbxLPS", "PostAbx_vs_PreAbxLPS")) {
    c(-1, 35)
  } else {
    c(-5, 5)
  }
  
  p <- ev_volcano(df, nm_pretty(nm), select_labels = sel, xlim_override = xlim_use)
  ggsave2(p, paste0(nm_pretty(nm), "_VP"), out_dir = volcano_dir, width = 5, height = 5)
}))

## 7.2 DEG counts barplot
count_data <- bind_rows(lapply(names(deg_list_sig), function(nm) {
  count_up_down(deg_list_sig[[nm]], nm)
}))

count_data <- count_data %>%
  mutate(
    Dataset_lbl = coalesce(label_map[Dataset], Dataset),
    Dataset_lbl = factor(Dataset_lbl, levels = rev(c(
      "Post-Abx vs Pre-Abx",
      "Pre-Abx + LPS vs Pre-Abx",
      "Post-Abx + LPS vs Pre-Abx + LPS",
      "Post-Abx + LPS vs Post-Abx",
      "Post-Abx vs Pre-Abx + LPS"
    ))),
    Status = factor(Status, levels = c("Down", "Up"))
  )

p_deg <- ggplot(count_data, aes(x = Dataset_lbl, y = Count, fill = Status)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  geom_text(
    aes(label = Count),
    position = position_dodge(width = 0.7),
    vjust = 0.5, hjust = -0.1,
    size = 1.5,
    family = "Arial"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c(Down = "#327ebaff", Up = "#da6868ff")) +
  labs(title = "Differentially Expressed Genes", x = NULL, y = "Count", fill = NULL) +
  theme_minimal(base_size = 7) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 7, family = "Arial"),
    axis.title = element_text(size = 7, family = "Arial"),
    axis.text = element_text(size = 6, family = "Arial"),
    legend.text = element_text(size = 6, family = "Arial"),
    legend.position = "bottom",
    plot.margin = margin(5.5, 30, 5.5, 5.5)
  )

ggsave2(
  p_deg,
  "DEG_counts_barplot",
  width = 7.5,
  height = 4.5,
  out_dir = file.path(deg_overview_dir, "Count-Barplot")
)
ggsave(file.path(deg_overview_dir, "DEG.svg"), plot = p_deg, width = 8, height = 5, units = "cm")

## 7.3 Venn and overlap analysis
sig_sets <- lapply(names(deg_list_sig), function(nm) {
  df <- deg_list_sig[[nm]]
  .check_cols(df, c("gene_id"))
  unique(as.character(df$gene_id))
})
names(sig_sets) <- unname(label_map[names(deg_list_sig)])

venn_dir <- file.path(deg_overview_dir, "Venn-Diagram")
dir.create(venn_dir, showWarnings = FALSE, recursive = TRUE)

p_venn <- NULL

if (requireNamespace("ggvenn", quietly = TRUE)) {
  p_venn <- ggvenn::ggvenn(
    sig_sets,
    fill_color = c("#FDBF6F", "#C98895", "#A6CEE3", "#CAB2D6", "#B2DF8A"),
    stroke_size = 0.4,
    set_name_size = 2.3,
    text_size = 1.8,
    show_stats = "cp"
  ) +
    ggtitle("Venn: significant DEGs (padj<0.05 & |log2FC|≥1.5)") +
    theme(
      text = element_text(family = "Arial", size = 7),
      plot.title = element_text(family = "Arial", size = 7),
      legend.text = element_text(family = "Arial"),
      axis.text = element_text(family = "Arial")
    )
} else if (requireNamespace("VennDiagram", quietly = TRUE)) {
  png(file.path(venn_dir, "Venn_sig.png"), width = 2000, height = 1800, res = 220)
  grid::grid.newpage()
  vd <- VennDiagram::venn.diagram(
    x = sig_sets,
    filename = NULL,
    fill = c("#FDBF6F", "#C98895", "#A6CEE3", "#CAB2D6", "#B2DF8A"),
    alpha = 0.5,
    cat.cex = 1.2,
    cex = 1.2,
    main = "Venn: significant DEGs (padj<0.05 & |log2FC|≥1.5)"
  )
  grid::grid.draw(vd)
  dev.off()
}

if (inherits(p_venn, "gg")) {
  ggsave(file.path(venn_dir, "Venn_sig.svg"), p_venn, width = 7, height = 7, units = "cm")
  ggsave(file.path(venn_dir, "Venn_sig.tiff"), p_venn, width = 12, height = 6, units = "cm",
         dpi = 300, compression = "lzw")
}

up_sets <- lapply(deg_list_sig, function(df) {
  df %>% filter(log2FoldChange >= lfc_thr, padj < alpha) %>% pull(gene_id) %>% unique()
})
down_sets <- lapply(deg_list_sig, function(df) {
  df %>% filter(log2FoldChange <= -lfc_thr, padj < alpha) %>% pull(gene_id) %>% unique()
})

if (requireNamespace("ggvenn", quietly = TRUE)) {
  ggsave2(
    ggvenn::ggvenn(
      up_sets,
      fill_color = c("#BFD3C1", "#C7CEEA", "#FFDAC1", "#E2F0CB", "#EFC3CA"),
      stroke_size = 0.5,
      set_name_size = 4
    ) + ggtitle("Venn: Up-regulated DEGs"),
    "Venn_up",
    width = 6.5,
    height = 6,
    out_dir = venn_dir
  )
  
  ggsave2(
    ggvenn::ggvenn(
      down_sets,
      fill_color = c("#BFD3C1", "#C7CEEA", "#FFDAC1", "#E2F0CB", "#EFC3CA"),
      stroke_size = 0.5,
      set_name_size = 4
    ) + ggtitle("Venn: Down-regulated DEGs"),
    "Venn_down",
    width = 6.5,
    height = 6,
    out_dir = venn_dir
  )
}

all_sheets <- pairwise_sheets(sig_sets)
writexl::write_xlsx(all_sheets, path = file.path(venn_dir, "shared_pairwise_all.xlsx"))

if (exists("up_sets")) {
  up_sheets <- pairwise_sheets(up_sets)
  writexl::write_xlsx(up_sheets, path = file.path(venn_dir, "shared_pairwise_up.xlsx"))
}
if (exists("down_sets")) {
  down_sheets <- pairwise_sheets(down_sets)
  writexl::write_xlsx(down_sheets, path = file.path(venn_dir, "shared_pairwise_down.xlsx"))
}

write_exclusive_xlsx(sig_sets, file.path(venn_dir, "exclusive_per_contrast_all.xlsx"))
if (exists("up_sets")) {
  write_exclusive_xlsx(up_sets, file.path(venn_dir, "exclusive_per_contrast_up.xlsx"))
}
if (exists("down_sets")) {
  write_exclusive_xlsx(down_sets, file.path(venn_dir, "exclusive_per_contrast_down.xlsx"))
}

############################################################
# 8. Shared-gene GO analysis
############################################################

go_shared_dir <- "GO_shared"
dir.create(go_shared_dir, showWarnings = FALSE, recursive = TRUE)

run_go_bp <- function(entrez, title_tag, file_tag, out_dir = go_shared_dir) {
  if (!length(entrez)) return(invisible(NULL))
  go <- tryCatch(
    enrichGO(
      gene = entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
      ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(go) || nrow(as.data.frame(go)) == 0) {
    writexl::write_xlsx(
      list(Results = data.frame(Message = paste("No GO:BP terms for", title_tag))),
      file.path(out_dir, paste0("GO_", file_tag, ".xlsx"))
    )
    return(invisible(NULL))
  }
  
  dp <- dotplot(go, showCategory = 10) + ggtitle(paste("GO:BP —", title_tag))
  ggsave2(dp, paste0("GO_", file_tag), out_dir = out_dir)
  writexl::write_xlsx(as.data.frame(go), file.path(out_dir, paste0("GO_", file_tag, ".xlsx")))
  invisible(go)
}

ct_a <- "PostAbx_vs_PreAbx"
ct_b <- "PreAbxLPS_vs_PreAbx"

genes_a <- unique(as.character(deg_list_sig[[ct_a]]$gene_id))
genes_b <- unique(as.character(deg_list_sig[[ct_b]]$gene_id))
shared_all <- sort(intersect(genes_a, genes_b))

entrez_all <- sym_to_entrez(shared_all)
run_go_bp(
  entrez_all,
  title_tag = "Shared DEGs: Post-Abx vs Pre-Abx ∩ Pre-Abx+LPS vs Pre-Abx",
  file_tag = "shared_all"
)

readr::write_tsv(
  tibble::tibble(gene_symbol = shared_all),
  file.path(go_shared_dir, "shared_gene_symbols.tsv")
)

############################################################
# 9. Functional enrichment analysis
############################################################

ora_results <- list()

for (ct in names(DEG_res_list)) {
  cat("Processing:", ct, "\n")
  df <- DEG_res_list[[ct]]
  
  ct_clean <- sanitize_name(ct, 80)
  ct_dir <- file.path(enrich_dir, ct_clean)
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
  
  sig_subsets <- list(
    All  = df |> filter(padj < alpha & abs(log2FoldChange) >= lfc_thr) |> pull(gene_id) |> unique(),
    Up   = df |> filter(padj < alpha &  log2FoldChange >=  lfc_thr)    |> pull(gene_id) |> unique(),
    Down = df |> filter(padj < alpha &  log2FoldChange <= -lfc_thr)    |> pull(gene_id) |> unique()
  )
  
  for (nm in names(sig_subsets)) {
    syms <- sig_subsets[[nm]]
    nm_clean <- sanitize_name(nm, 40)
    sub_dir <- file.path(ct_dir, nm_clean)
    dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)
    
    if (length(syms) == 0) next
    
    entrez <- sym_to_entrez(syms)
    if (length(entrez) == 0) next
    
    ora <- run_ora(entrez)
    
    if (is.null(ora_results[[ct]])) ora_results[[ct]] <- list()
    ora_results[[ct]][[nm]] <- list(
      symbols = syms,
      entrez  = entrez,
      GO      = ora$GO,
      KEGG    = ora$KEGG,
      REACT   = ora$REACT
    )
    
    write_or_placeholder(ora$GO,    file.path(sub_dir, paste0("GO_",    ct_clean, "_", nm_clean, ".xlsx")),    paste(ct, nm, "GO"))
    write_or_placeholder(ora$KEGG,  file.path(sub_dir, paste0("KEGG_",  ct_clean, "_", nm_clean, ".xlsx")),    paste(ct, nm, "KEGG"))
    write_or_placeholder(ora$REACT, file.path(sub_dir, paste0("React_", ct_clean, "_", nm_clean, ".xlsx")),    paste(ct, nm, "Reactome"))
  }
}

saveRDS(ora_results, file.path(enrich_dir, "ORA_results.rds"))

go_post_vs_pre <- ora_results[["PostAbx_vs_PreAbx"]][["All"]][["GO"]]@result

terms_exact <- c(
  "cellular response to lipopolysaccharide",
  "response to tumor necrosis factor",
  "chemotaxis",
  "granulocyte migration",
  "regulation of inflammatory response"
)

p_go_post_vs_pre <- go_dotplot_df(
  df = go_post_vs_pre,
  include_desc = terms_exact,
  regex = FALSE,
  title = "DEGs - Post-Abx vs Pre-Abx (GO:BP)",
  out_file = file.path(enrich_dir, "Plots", "GO", "PostAbx_vs_PreAbx")
)

## 9.1 compareCluster
sets_all <- lapply(DEG_res_list, \(d) d |> filter(padj < alpha & abs(log2FoldChange) >= lfc_thr) |> pull(gene_id) |> unique())
sets_up  <- lapply(DEG_res_list, \(d) d |> filter(padj < alpha &  log2FoldChange >=  lfc_thr)    |> pull(gene_id) |> unique())
sets_dn  <- lapply(DEG_res_list, \(d) d |> filter(padj < alpha &  log2FoldChange <= -lfc_thr)    |> pull(gene_id) |> unique())

names(sets_all) <- names(DEG_res_list)
names(sets_up)  <- names(DEG_res_list)
names(sets_dn)  <- names(DEG_res_list)

cc_dir <- file.path(enrich_dir, "Plots", "compareCluster")
dir.create(cc_dir, showWarnings = FALSE, recursive = TRUE)

cc_all  <- run_cc(sets_all, "All", cc_dir)
cc_up   <- run_cc(sets_up,  "Up",  cc_dir)
cc_down <- run_cc(sets_dn,  "Down", cc_dir)

## 9.2 Network analysis
results_network <- list()

for (nm in names(deg_list_sig)) {
  message("\n### Processing contrast:", nm, "###")
  df <- deg_list_sig[[nm]]
  
  results_network[[nm]] <- list(
    All  = perform_network_analysis(df, nm, "All"),
    Up   = perform_network_analysis(filter(df, padj < 0.05, log2FoldChange > 0), nm, "Up"),
    Down = perform_network_analysis(filter(df, padj < 0.05, log2FoldChange < 0), nm, "Down")
  )
}


############################################################
# 10. TF and cofactor analysis
############################################################

tf_summary_dir  <- file.path(tf_dir, "TF_DEG_Summaries")
cof_summary_dir <- file.path(tf_dir, "CoF_DEG_Summaries")

TF_results <- extract_deg_class(
  datasets_filt = deg_list_sig,
  class_file = "Homo_sapiens_TF.txt",
  out_dir = tf_summary_dir,
  class_label = "TF"
)

CoF_results <- extract_deg_class(
  datasets_filt = deg_list_sig,
  class_file = "Homo_sapiens_CoF.txt",
  out_dir = cof_summary_dir,
  class_label = "CoF"
)

for (contrast in names(TF_results)) {
  plot_deg_class(TF_results[[contrast]]$All, contrast, "TF", out_dir = file.path(tf_summary_dir, "Plot"))
}

for (contrast in names(CoF_results)) {
  plot_deg_class(CoF_results[[contrast]]$All, contrast, "CoF", out_dir = file.path(cof_summary_dir, "Plot"))
}

############################################################
# 11. Export summary files
############################################################

PostAbx_vs_PreAbx <- deg_list_sig[["PostAbx_vs_PreAbx"]]
save(PostAbx_vs_PreAbx, file = file.path(deg_tables_dir, "DEG_PostAbx_vs_PreAbx_sig_FC1.5.RData"))

save_list_to_excel(DEG_res_list, output_file = file.path(deg_tables_dir, "DEGs_All.xlsx"))
save_list_to_excel(deg_list_sig, output_file = file.path(deg_tables_dir, "DEGs_Filtered.xlsx"))


############################################################
# 12. Save full workspace
############################################################

save.image(file = file.path(out_dir, "DEG-plots.RData"))



############################################################
# End of script
############################################################