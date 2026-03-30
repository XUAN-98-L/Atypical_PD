# Pathway enrichment figure: terms + significance + LogFC + function (one per contrast)
# Reads all GO_BP_ssGSEA_*_all_terms.csv from Results/01_Top10_6comparisons

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)
library(fmsb)
library(ComplexHeatmap)
library(circlize)
library(grid)

my_palette = c("#66C2A5","#E78AC3", "#ffc773",
            "#1B9E77","#7570B3","#FEE090",
            "#A6D854","#8DA0CB","#FC8D62",
            "#8FBC94","#b0a4e3","#ffa631",
            "#0aa344","#e4c6d0","#ffa400",
            "#519a73","#4b5cc4","#eedeb0",
            "#549688","#ffb3a7","#b35c44",
            "#7fecad","#a1afc9","#a78e44",
            "#519a73","#2e4e7e","#955539")

# Paths
input_dir <- "Results/01_Top10_6comparisons"
output_dir <- "Results/03_compare_control"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Find and read all *_all_terms.csv files
all_terms_files <- list.files(input_dir, pattern = "_all_terms\\.csv$", full.names = TRUE)
dat_list <- lapply(all_terms_files, function(f) {
  contrast <- basename(f) %>%
    gsub("GO_BP_ssGSEA_", "", .) %>%
    gsub("_all_terms\\.csv$", "", .)
  read_csv(f, show_col_types = FALSE) %>% mutate(contrast = contrast)
})
dat_full <- bind_rows(dat_list)

# Top10 list (term, contrast) for filtering; keep only vs Ctr comparisons
top10_path <- file.path(input_dir, "GO_BP_ssGSEA_top10_terms_for_heatmap.csv")
top10_raw <- read_csv(top10_path, show_col_types = FALSE) %>% select(term, contrast)

top10_list <- top10_raw %>%
  filter(contrast %in% c("NPS_vs_Ctr", "Tauopathy_vs_Ctr", "PD_vs_Ctr"))

# Top10 union across these five contrasts (see table(top10_list_all_contrasts$contrast))
all_contrasts_heatmap <- c(
  "NPS_vs_Ctr", "NPS_vs_PD", "PD_vs_Ctr", "Tauopathy_vs_Ctr", "Tauopathy_vs_PD"
)
top10_list_all_contrasts <- top10_raw %>%
  filter(contrast %in% all_contrasts_heatmap)
terms_all_contrasts <- unique(top10_list_all_contrasts$term)

# Format term label: remove GOBP_ prefix, replace _ with space, title case
format_term <- function(x) {
  x %>%
    gsub("^GOBP_", "", .) %>%
    gsub("_", " ", .) %>%
    tolower() %>%
    tools::toTitleCase()
}

# Optional: provide a manual mapping table with columns: term, Function
# - term should match the GO term IDs (e.g., "GOBP_...") used in the CSVs
# - Function is your category label (any string)
function_map_path <- "Data/Mannual_annotation.xlsx"
library(readxl)
function_map <- read_excel(function_map_path) %>%
  rename(Function = any_of(c("Function", "function"))) %>%
  select(term, Function) %>%
  distinct()

get_function <- function(term) {
  hit <- function_map$Function[match(term, function_map$term)]
  if (is.na(hit) || !nzchar(hit)) {
    stop(
      "Missing Function mapping for term: ", term,
      "\nPlease add it to: ", function_map_path
    )
  }
  hit
}

# Build wide top10_terms: one row per term, columns term, group, {contrast}_logFC, {contrast}_adj.P.Val
vs_ctr <- c("NPS_vs_Ctr", "PD_vs_Ctr", "Tauopathy_vs_Ctr")
dat_vs_ctr <- dat_full %>%
  filter(term %in% unique(top10_list$term), contrast %in% vs_ctr)
logfc_wide <- dat_vs_ctr %>%
  select(term, contrast, logFC) %>%
  pivot_wider(names_from = contrast, values_from = logFC, names_glue = "{contrast}_logFC")
adjp_wide <- dat_vs_ctr %>%
  select(term, contrast, adj.P.Val) %>%
  pivot_wider(names_from = contrast, values_from = adj.P.Val, names_glue = "{contrast}_adj.P.Val")
top10_terms <- logfc_wide %>%
  full_join(adjp_wide, by = "term") #%>%
#   mutate(group = sapply(term, assign_function)) %>%
#   relocate(group, .after = term)

# Add manual/derived function group and keep contrast from top10 list
top10_terms <- top10_terms %>%
  mutate(group = sapply(term, get_function)) %>%
  left_join(top10_list, by = "term")

# Export wide top10_terms table
write.csv(top10_terms, file.path(output_dir, "GO_BP_ssGSEA_top10_terms_wide.csv"), row.names = FALSE, quote = FALSE)

# Wide table: all five contrasts (union of top10 terms per contrast), same layout as above
dat_all_contrasts_tbl <- dat_full %>%
  filter(term %in% terms_all_contrasts, contrast %in% all_contrasts_heatmap)
logfc_wide_all <- dat_all_contrasts_tbl %>%
  select(term, contrast, logFC) %>%
  pivot_wider(names_from = contrast, values_from = logFC, names_glue = "{contrast}_logFC")
adjp_wide_all <- dat_all_contrasts_tbl %>%
  select(term, contrast, adj.P.Val) %>%
  pivot_wider(names_from = contrast, values_from = adj.P.Val, names_glue = "{contrast}_adj.P.Val")
top10_contrast_tags <- top10_list_all_contrasts %>%
  distinct(term, contrast) %>%
  group_by(term) %>%
  summarise(contrast = paste(contrast, collapse = ";"), .groups = "drop")

top10_terms_all_contrasts <- logfc_wide_all %>%
  full_join(adjp_wide_all, by = "term") %>%
  mutate(group = sapply(term, get_function)) %>%
  left_join(top10_contrast_tags, by = "term")

write.csv(
  top10_terms_all_contrasts,
  file.path(output_dir, "GO_BP_ssGSEA_top10_terms_wide_all_contrasts.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Contrast label for plot (e.g. "NPS vs PD")
contrast_label <- function(c) {
  c %>%
    gsub("_vs_", " vs ", .) %>%
    gsub("_", " ", .)
}

# Function colors (use your custom palette)
# Keep Function order as it appears in your manual annotation table.
fun_levels <- unique(function_map$Function)
if (length(my_palette) < length(fun_levels)) {
  stop(
    "my_palette has ", length(my_palette), " colors, but you have ",
    length(fun_levels), " Function levels. Please add more colors to my_palette."
  )
}
fun_cols <- my_palette[seq_along(fun_levels)]
fun_colors <- setNames(fun_cols, fun_levels)

# Combined heatmap for 3 comparisons:
# - color = logFC
# - stars = P.Value thresholds
vs_ctr <- c("NPS_vs_Ctr", "PD_vs_Ctr", "Tauopathy_vs_Ctr")
plot_df <- dat_full %>%
  filter(term %in% unique(top10_list$term), contrast %in% vs_ctr) %>%
  mutate(
    term_label = format_term(term),
    Function = sapply(term, get_function),
    sig = case_when(
      P.Value < 0.0001 ~ "****",
      P.Value < 0.001  ~ "***",
      P.Value < 0.01   ~ "**",
      P.Value < 0.05   ~ "*",
      TRUE ~ ""
    )
  )

logfc_mat <- plot_df %>%
  select(term_label, contrast, logFC) %>%
  pivot_wider(names_from = contrast, values_from = logFC) %>%
  column_to_rownames("term_label") %>%
  as.matrix()

sig_mat <- plot_df %>%
  select(term_label, contrast, sig) %>%
  pivot_wider(names_from = contrast, values_from = sig) %>%
  column_to_rownames("term_label") %>%
  as.matrix()

# Prettify column names for plotting
colnames(logfc_mat) <- gsub("_", " ", colnames(logfc_mat))
colnames(sig_mat) <- gsub("_", " ", colnames(sig_mat))

# Copy for clustered heatmap (no Function order / split)
logfc_mat_for_cluster <- logfc_mat
sig_mat_for_cluster <- sig_mat

fun_vec <- plot_df %>%
  distinct(term_label, Function) %>%
  tibble::deframe()

# Order rows like the top10_list appearance (within each contrast) as much as possible
row_order <- top10_list %>%
  filter(contrast %in% vs_ctr) %>%
  mutate(term_label = format_term(term)) %>%
  pull(term_label) %>%
  unique()
row_order <- row_order[row_order %in% rownames(logfc_mat)]
logfc_mat <- logfc_mat[row_order, , drop = FALSE]
sig_mat <- sig_mat[row_order, , drop = FALSE]
fun_vec <- fun_vec[row_order]

# Keep Function groups together, but hide the slice titles text
fun_split <- factor(fun_vec, levels = unique(fun_vec))
fun_split_titles <- rep("", length(levels(fun_split)))

# Color scale centered at 0 (use full term set so scales match both heatmaps)
fc_range <- range(logfc_mat_for_cluster, na.rm = TRUE)
fc_max <- max(abs(fc_range))
col_logfc <- circlize::colorRamp2(c(-fc_max, 0, fc_max), c("navy", "white", "firebrick3"))

ha_left <- ComplexHeatmap::rowAnnotation(
  Function = fun_vec,
  col = list(Function = fun_colors),
  show_annotation_name = FALSE,
  width = unit(4, "mm")
)

ht <- ComplexHeatmap::Heatmap(
  logfc_mat,
  name = "logFC",
  col = col_logfc,
  row_split = fun_split,
  row_title = fun_split_titles,
  row_title_gp = gpar(fontsize = 0),
  row_title_rot = 0,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 10),
  cell_fun = function(j, i, x, y, width, height, fill) {
    lab <- sig_mat[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      grid.text(lab, x, y, gp = gpar(fontsize = 8, fontface = "bold", col = "black"))
    }
  },
  heatmap_legend_param = list(
    title = "logFC",
    direction = "vertical"
  )
)

pdf(file.path(output_dir, "GO_BP_ssGSEA_top10_3comparisons_logFC_with_stars.pdf"),
    width = 7, height = max(4, nrow(logfc_mat) * 0.25))
ComplexHeatmap::draw(
  ha_left + ht,
  annotation_legend_side = "right",
  heatmap_legend_side = "right",
  padding = unit(c(2, 6, 2, 2), "mm")
)
dev.off()

# Same heatmap without Function annotation; rows clustered (no Function split / ordering).
# Use cluster_rows = dendrogram (Heatmap has no row_dend argument; partial match breaks).
hc_row <- stats::hclust(stats::dist(logfc_mat_for_cluster))

ht_cluster <- ComplexHeatmap::Heatmap(
  logfc_mat_for_cluster,
  name = "logFC",
  col = col_logfc,
  cluster_rows = stats::as.dendrogram(hc_row),
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 10),
  cell_fun = function(j, i, x, y, width, height, fill) {
    lab <- sig_mat_for_cluster[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      grid.text(lab, x, y, gp = gpar(fontsize = 8, fontface = "bold", col = "black"))
    }
  },
  heatmap_legend_param = list(
    title = "logFC",
    direction = "vertical"
  )
)

pdf(file.path(output_dir, "GO_BP_ssGSEA_top10_3comparisons_logFC_with_stars_rowCluster_noFunction.pdf"),
    width = 7, height = max(4, nrow(logfc_mat_for_cluster) * 0.25))
ComplexHeatmap::draw(
  ht_cluster,
  heatmap_legend_side = "right",
  padding = unit(c(2, 6, 2, 2), "mm")
)
dev.off()

# Heatmap: union of top10 terms across NPS_vs_Ctr, NPS_vs_PD, PD_vs_Ctr,
# Tauopathy_vs_Ctr, Tauopathy_vs_PD — no Function, row dendrogram only
plot_df_all <- dat_full %>%
  filter(term %in% terms_all_contrasts, contrast %in% all_contrasts_heatmap) %>%
  mutate(
    term_label = format_term(term),
    sig = case_when(
      P.Value < 0.0001 ~ "****",
      P.Value < 0.001  ~ "***",
      P.Value < 0.01   ~ "**",
      P.Value < 0.05   ~ "*",
      TRUE ~ ""
    )
  )

logfc_mat_all <- plot_df_all %>%
  select(term_label, contrast, logFC) %>%
  pivot_wider(names_from = contrast, values_from = logFC) %>%
  column_to_rownames("term_label") %>%
  as.matrix()
logfc_mat_all <- logfc_mat_all[, all_contrasts_heatmap, drop = FALSE]

sig_mat_all <- plot_df_all %>%
  select(term_label, contrast, sig) %>%
  pivot_wider(names_from = contrast, values_from = sig) %>%
  column_to_rownames("term_label") %>%
  as.matrix()
sig_mat_all <- sig_mat_all[, all_contrasts_heatmap, drop = FALSE]

colnames(logfc_mat_all) <- gsub("_", " ", colnames(logfc_mat_all))
colnames(sig_mat_all) <- gsub("_", " ", colnames(sig_mat_all))

if (any(is.na(logfc_mat_all))) {
  stop("NA in logFC matrix for 5-contrast heatmap; check all_terms CSVs.")
}

fc_range_all <- range(logfc_mat_all, na.rm = TRUE)
fc_max_all <- max(abs(fc_range_all))
col_logfc_all <- circlize::colorRamp2(
  c(-fc_max_all, 0, fc_max_all),
  c("navy", "white", "firebrick3")
)

hc_row_all <- stats::hclust(stats::dist(logfc_mat_all))

ht_all <- ComplexHeatmap::Heatmap(
  logfc_mat_all,
  name = "logFC",
  col = col_logfc_all,
  cluster_rows = stats::as.dendrogram(hc_row_all),
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  cell_fun = function(j, i, x, y, width, height, fill) {
    lab <- sig_mat_all[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      grid.text(lab, x, y, gp = gpar(fontsize = 7, fontface = "bold", col = "black"))
    }
  },
  heatmap_legend_param = list(
    title = "logFC",
    direction = "vertical"
  )
)

pdf(
  file.path(output_dir, "GO_BP_ssGSEA_top10_5comparisons_logFC_with_stars_rowCluster_noFunction.pdf"),
  width = 9,
  height = max(5, nrow(logfc_mat_all) * 0.22)
)
ComplexHeatmap::draw(
  ht_all,
  heatmap_legend_side = "right",
  padding = unit(c(2, 6, 2, 2), "mm")
)
dev.off()

# Rotated row dendrogram: put two cutree clusters at the bottom of the heatmap and
# swap their order (purple vs red). Pick k and cluster IDs after coloring the tree once:
#   d0 <- as.dendrogram(hc_row_all)
#   plot(dendextend::color_branches(d0, k = 5))
suppressPackageStartupMessages(library(dendextend))

rotate_row_k <- 6L
Immune1 <- 4L
Protein_metabolism <- 2L
phophorylation <- 1L
other_cellular_processes <- 3L
Immune2 <- 5L
Chemotaxis <- 6L

d_row <- as.dendrogram(hc_row_all)
hc_row_all <- stats::hclust(stats::dist(logfc_mat_all))
grp <- stats::cutree(hc_row_all, k = rotate_row_k)
leaf_in_dend_order <- rownames(logfc_mat_all)[hc_row_all$order]

other <- leaf_in_dend_order[
  !grp[leaf_in_dend_order] %in% c(phophorylation, other_cellular_processes)
]

blk_Immune1 <- leaf_in_dend_order[grp[leaf_in_dend_order] == Immune1]
blk_Immune2 <- leaf_in_dend_order[grp[leaf_in_dend_order] == Immune2]
blk_Protein_metabolism <- leaf_in_dend_order[grp[leaf_in_dend_order] == Protein_metabolism]
blk_phophorylation <- leaf_in_dend_order[grp[leaf_in_dend_order] == phophorylation]
blk_other_cellular_processes <- leaf_in_dend_order[grp[leaf_in_dend_order] == other_cellular_processes]
blk_Chemotaxis <- leaf_in_dend_order[grp[leaf_in_dend_order] == Chemotaxis]

new_leaf_order <- c(blk_Immune1, blk_Immune2, blk_Chemotaxis, blk_Protein_metabolism, blk_phophorylation, blk_other_cellular_processes,other)

d_row_rot <- dendextend::rotate(d_row, order = new_leaf_order)

ht_all_rot <- ComplexHeatmap::Heatmap(
  logfc_mat_all,
  name = "logFC",
  col = col_logfc_all,
  cluster_rows = d_row_rot,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  cell_fun = function(j, i, x, y, width, height, fill) {
    lab <- sig_mat_all[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      grid.text(lab, x, y, gp = gpar(fontsize = 7, fontface = "bold", col = "black"))
    }
  },
  heatmap_legend_param = list(
    title = "logFC",
    direction = "vertical"
  )
)

pdf(
  file.path(output_dir, "GO_BP_ssGSEA_top10_5comparisons_logFC_with_stars_rowCluster_noFunction_rotate.pdf"),
  width = 9,
  height = max(5, nrow(logfc_mat_all) * 0.22)
)
ComplexHeatmap::draw(
  ht_all_rot,
  heatmap_legend_side = "right",
  padding = unit(c(2, 6, 2, 2), "mm")
)
dev.off()

# Spider / radar plot (mean abs(logFC) per Function; normalized to [0, 1])
# Make it comparable and not "weird":
# - fix Function order
# - complete missing Function for each comparison with 0
# - draw on numeric angle axis and close explicitly
radar_df <- dat_full %>%
  filter(term %in% unique(top10_list$term), contrast %in% vs_ctr) %>%
  mutate(
    Function = sapply(term, get_function),
    contrast_label = gsub("_", " ", contrast),
    score_raw = abs(logFC)
  ) %>%
  group_by(Function, contrast_label) %>%
  summarise(score_raw = mean(score_raw, na.rm = TRUE), .groups = "drop") %>%
  mutate(Function = factor(Function, levels = fun_levels)) %>%
  tidyr::complete(Function, contrast_label, fill = list(score_raw = 0))

max_abs_logfc <- max(radar_df$score_raw, na.rm = TRUE)
if (!is.finite(max_abs_logfc) || max_abs_logfc <= 0) {
  stop("max abs(logFC) is not positive; cannot normalize radar plot.")
}
radar_df <- radar_df %>% mutate(score = score_raw / max_abs_logfc)

radar_df <- radar_df %>% mutate(angle_id = as.integer(Function))

radar_df_closed <- radar_df %>%
  arrange(contrast_label, angle_id) %>%
  group_by(contrast_label) %>%
  group_modify(~ bind_rows(.x, dplyr::slice(.x, 1) %>% mutate(angle_id = max(.x$angle_id) + 1L))) %>%
  ungroup()

contrast_levels <- unique(radar_df$contrast_label)
# Colors follow Script/01_Top10_6comparisons.R (112-115)
comparison_colors <- c(
  "NPS vs Ctr" = "#8DA0CB",
  "PD vs Ctr" = "#A6D854",
  "Tauopathy vs Ctr" = "#E78AC3"
)
contrast_cols <- unname(comparison_colors[contrast_levels])
if (any(is.na(contrast_cols))) {
  stop("Missing color mapping for: ", paste(contrast_levels[is.na(contrast_cols)], collapse = ", "))
}

# Build fmsb radarchart input:
# first two rows are max/min, following rows are each comparison
radar_wide <- radar_df %>%
  select(Function, contrast_label, score) %>%
  pivot_wider(names_from = Function, values_from = score) %>%
  arrange(factor(contrast_label, levels = contrast_levels))

radar_mat <- as.data.frame(radar_wide %>% select(-contrast_label))
radar_mat <- rbind(rep(1, ncol(radar_mat)), rep(0, ncol(radar_mat)), radar_mat)
radar_mat <- as.data.frame(radar_mat)
colnames(radar_mat) <- colnames(radar_wide %>% select(-contrast_label))

pdf(file.path(output_dir, "GO_BP_spider_meanAbsLogFC_byFunction_norm.pdf"),
    width = 10, height = 6)
fmsb::radarchart(
  radar_mat,
  axistype = 1,
  pcol = contrast_cols,
  pfcol = grDevices::adjustcolor(contrast_cols, alpha.f = 0.25),
  plwd = 2,
  plty = 1,
  cglcol = "grey85",
  cglty = 1,
  axislabcol = "grey40",
  caxislabels = c("0", "0.25", "0.5", "0.75", "1.0"),
  vlcex = 0.9
)

# Caption/annotation for interpretation
title(
  main = "Mean abs(logFC) per Function (top10 terms), normalized to [0, 1]"
)
legend(
  "topright",
  legend = contrast_levels,
  bty = "n",
  pch = 16,
  col = contrast_cols,
  pt.cex = 1.2,
  cex = 1
)
dev.off()

message("Figures saved to ", output_dir)
