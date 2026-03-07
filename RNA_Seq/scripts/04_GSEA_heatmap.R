############################################################
# Script: 04_GSEA_heatmap.R
# Purpose: Generate heatmap visualizations for selected
#          genes from GSEA-enriched signatures.
# Input: 6_GSEA/, 1_DEG_tables/
# Output: 6_GSEA/Heatmap_inflammation.svg
# Author: Bruna Toledo
############################################################

############################################################
# 1. Load libraries
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

############################################################
# 2. Define input/output paths
############################################################

gsea_dir      <- "6_GSEA"
deg_tables_dir <- "1_DEG_tables"

hallmark_file <- file.path(gsea_dir, "GENE_SYMBOLS_GSEA_INFLAMMATORY.txt")
output_file   <- file.path(gsea_dir, "Heatmap_inflammation.svg")

############################################################
# 3. Check required objects
############################################################

required_objects <- c("deg_list_sig", "norm_counts_anno", "metadata")

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "The following objects are missing from the environment: ",
    paste(missing_objects, collapse = ", "),
    ". Please load the DEG/count workspace before running this script."
  )
}

############################################################
# 4. Load hallmark gene set
############################################################

hall_inf <- read.delim(hallmark_file, stringsAsFactors = FALSE)
hall_inf <- hall_inf$GENE_SYMBOLS

############################################################
# 5. Select genes of interest
############################################################

PostAbx_vs_PreAbx_sig <- deg_list_sig[["PostAbx_vs_PreAbx"]]

deg_hall_inf <- PostAbx_vs_PreAbx_sig %>%
  filter(gene_id %in% hall_inf) %>%
  pull(gene_id)

############################################################
# 6. Extract normalized expression matrix
############################################################

norm_filt <- norm_counts_anno %>%
  filter(gene_name %in% deg_hall_inf)

norm_filt_num <- norm_filt
rownames(norm_filt_num) <- norm_filt_num$gene_name
norm_filt_num <- norm_filt_num[, c(3:14)]
norm_filt_num <- as.matrix(norm_filt_num)

############################################################
# 7. Subset samples for heatmap
############################################################

norm_filt_num_post <- norm_filt_num[, c(1, 3, 5, 7, 9, 11)]

anno_col_post <- data.frame(
  Group   = metadata$condition,
  Subject = metadata$subject
)

anno_col_post <- as.data.frame(anno_col_post[c(1, 3, 5, 7, 9, 11), ])
colnames(anno_col_post) <- c("Condition", "Subject")

anno_col_post$Condition <- factor(
  anno_col_post$Condition,
  levels = c("Pre-Abx", "Post-Abx")
)

anno_col_post$Group_plot <- factor(
  ifelse(anno_col_post$Condition == "Pre-Abx", "A Pre-Abx", "B Post-Abx"),
  levels = c("A Pre-Abx", "B Post-Abx")
)

############################################################
# 8. Build heatmap annotations
############################################################

anno_col <- HeatmapAnnotation(
  Group = anno_col_post$Group_plot,
  col = list(
    Group = c(
      "A Pre-Abx"  = "#804D8066",
      "B Post-Abx" = "#804D80cc"
    )
  ),
  annotation_name_gp = gpar(fontfamily = "Arial", fontsize = 6),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Group = list(
      title = "",
      title_gp = gpar(fontfamily = "Arial", fontsize = 6),
      labels_gp = gpar(fontfamily = "Arial", fontsize = 6)
    )
  )
)

spl <- droplevels(anno_col_post$Group_plot)
n_slices <- nlevels(spl)

############################################################
# 9. Scale expression matrix
############################################################

mat_z <- t(apply(norm_filt_num_post, 1, scale))

colnames(mat_z) <- c("S8", "S8",
                     "S9", "S9",
                     "S11", "S11")

color_z <- colorRamp2(
  seq(min(mat_z), max(mat_z), length = 3),
  c("royalblue2", "#EEEEEE", "firebrick2"),
  space = "RGB"
)

############################################################
# 10. Generate heatmap
############################################################

h_z <- Heatmap(
  mat_z,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  col = color_z,
  column_labels = colnames(mat_z),
  name = "Z-score",
  top_annotation = anno_col,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 5.5, fontfamily = "Arial"),
  column_names_gp = gpar(fontsize = 6, fontfamily = "Arial"),
  column_split = spl,
  column_title = rep("", n_slices),
  show_parent_dend_line = FALSE,
  show_column_names = FALSE,
  clustering_distance_rows = "manhattan",
  clustering_distance_columns = "manhattan",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    color_bar = "continuous",
    legend_width = unit(5, "cm"),
    title_position = "lefttop",
    title_gp = gpar(fontfamily = "Arial", fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontfamily = "Arial", fontsize = 6)
  )
)

############################################################
# 11. Draw heatmap
############################################################

hmap_z <- draw(
  h_z,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE
)

############################################################
# 12. Export figure
############################################################

svg(output_file, width = 2, height = 3.8)

draw(
  h_z,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE
)

dev.off()

############################################################
# 13. Save session information
############################################################

sink(paste0(out_prefix, "_sessionInfo.txt"))

print(sessionInfo())

sink()



############################################################
# End of script
############################################################