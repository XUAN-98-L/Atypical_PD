# Rscript Script/04_heatmap_vis_rotate.R -p 50 -k 6 --cluster-last A

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dendextend)
  library(GO.db)
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(optparse)
})

my_palette = c("#66C2A5","#E78AC3", "#ffc773",
            "#1B9E77","#7570B3","#FEE090",
            "#A6D854","#8DA0CB","#FC8D62",
            "#8FBC94","#b0a4e3","#ffa631",
            "#0aa344","#e4c6d0","#ffa400",
            "#519a73","#4b5cc4","#eedeb0",
            "#549688","#ffb3a7","#b35c44",
            "#7fecad","#a1afc9","#a78e44",
            "#519a73","#2e4e7e","#955539")

option_list <- list(
  optparse::make_option(
    c("-i", "--input-dir"),
    type = "character",
    default = "Results/01_Top10_6comparisons",
    dest = "input_dir",
    metavar = "PATH",
    help = "Directory of GO_BP_ssGSEA * _all_terms.csv and top10-for-heatmap CSV. Default: %default"
  ),
  optparse::make_option(
    c("-o", "--output-dir"),
    type = "character",
    default = "Results/04_heatmap_vis_rotate",
    dest = "output_dir",
    metavar = "PATH",
    help = "Output directory (PDF, exported tables). Default: %default"
  ),
  optparse::make_option(
    c("-a", "--row-annot-csv"),
    type = "character",
    default = "Data/GO_BP_ssGSEA_top10_5comparisons_rotate_cluster_IC_per_term.csv",
    dest = "row_annot_csv",
    metavar = "FILE",
    help = paste0(
      "Row-order + per-row Cluster + IC (term_label, term, GO_id, IC, Cluster). ",
      "Cluster defines Function/gaps; not cut from the dendrogram. Default: %default"
    )
  ),
  optparse::make_option(
    c("-p", "--percentage"),
    type = "double",
    default = 50,
    dest = "percentage",
    metavar = "NUM",
    help = paste0(
      "Minimum fraction of cluster GO terms (in percent) that must share a BP ancestor ",
      "for that ancestor to be listed. Example: 60 means >=60%%. Default: %default"
    )
  ),
  optparse::make_option(
    c("-k", "--rotate-row-k"),
    type = "integer",
    default = 6L,
    dest = "rotate_row_k",
    metavar = "INT",
    help = paste0(
      "k used by cutree() to build row blocks before dendrogram rotation. ",
      "Blocks are re-ordered by CSV Cluster order. Default: %default"
    )
  ),
  optparse::make_option(
    c("--cluster-last"),
    type = "character",
    default = "A",
    dest = "cluster_last",
    metavar = "CLUSTER_ID",
    help = paste0(
      "Optional: move this Cluster to the bottom in rotated row order. ",
      "Set empty string to disable. Default: %default"
    )
  )
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
input_dir <- opt$input_dir
output_dir <- opt$output_dir
row_annot_csv_path <- opt$row_annot_csv
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
common_ancestor_pct <- as.double(opt$percentage)
if (length(common_ancestor_pct) != 1L || is.na(common_ancestor_pct) ||
  common_ancestor_pct <= 0 || common_ancestor_pct > 100) {
  stop("--percentage must be a single number in (0, 100].")
}
rotate_row_k <- as.integer(opt$rotate_row_k)
if (length(rotate_row_k) != 1L || is.na(rotate_row_k) || rotate_row_k < 2L) {
  stop("--rotate-row-k must be an integer >= 2.")
}
cluster_last <- as.character(opt$cluster_last)
if (length(cluster_last) != 1L || is.na(cluster_last)) {
  stop("--cluster-last must be a single character string (or empty to disable).")
}
common_ancestor_frac <- common_ancestor_pct / 100
pct_fname_tag <- gsub(".", "p", sprintf("%.4g", common_ancestor_pct), fixed = TRUE)

format_term <- function(x) {
  x %>%
    gsub("^GOBP_", "", .) %>%
    gsub("_", " ", .) %>%
    tolower() %>%
    tools::toTitleCase()
}

`%||%` <- function(x, y) if (is.null(x)) y else x

if (!file.exists(row_annot_csv_path)) {
  stop("Missing row-ordered annotation CSV: ", row_annot_csv_path)
}
row_annot_tbl <- readr::read_csv(row_annot_csv_path, show_col_types = FALSE)
if (!all(c("term_label", "term", "GO_id", "IC") %in% colnames(row_annot_tbl))) {
  stop("Annotation CSV must include at least: term_label, term, GO_id, IC")
}
if (!"Cluster" %in% colnames(row_annot_tbl) && !"cluster" %in% colnames(row_annot_tbl)) {
  stop("Annotation CSV must include column Cluster (row grouping; not dendrogram cutree).")
}
row_annot_tbl <- row_annot_tbl %>% mutate(IC = suppressWarnings(as.numeric(IC)))
if ("Cluster" %in% colnames(row_annot_tbl)) {
  row_annot_tbl$cluster <- as.character(row_annot_tbl$Cluster)
} else {
  row_annot_tbl$cluster <- as.character(row_annot_tbl$cluster)
}
if (any(is.na(row_annot_tbl$cluster) | !nzchar(row_annot_tbl$cluster))) {
  stop("Column Cluster has missing or empty values in annotation CSV.")
}

# Build the same 5-comparison matrices as 03_compare_control.R
all_terms_files <- list.files(input_dir, pattern = "_all_terms\\.csv$", full.names = TRUE)
dat_list <- lapply(all_terms_files, function(f) {
  contrast <- basename(f) %>%
    gsub("GO_BP_ssGSEA_", "", .) %>%
    gsub("_all_terms\\.csv$", "", .)
  readr::read_csv(f, show_col_types = FALSE) %>% mutate(contrast = contrast)
})
dat_full <- bind_rows(dat_list)

top10_path <- file.path(input_dir, "GO_BP_ssGSEA_top10_terms_for_heatmap.csv")
top10_raw <- readr::read_csv(top10_path, show_col_types = FALSE) %>% dplyr::select(term, contrast)
all_contrasts_heatmap <- c("NPS_vs_Ctr", "NPS_vs_PD", "PD_vs_Ctr", "Tauopathy_vs_Ctr", "Tauopathy_vs_PD")
terms_all_contrasts <- top10_raw %>%
  filter(contrast %in% all_contrasts_heatmap) %>%
  pull(term) %>%
  unique()

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
  dplyr::select(term_label, contrast, logFC) %>%
  pivot_wider(names_from = contrast, values_from = logFC) %>%
  column_to_rownames("term_label") %>%
  as.matrix()
sig_mat_all <- plot_df_all %>%
  dplyr::select(term_label, contrast, sig) %>%
  pivot_wider(names_from = contrast, values_from = sig) %>%
  column_to_rownames("term_label") %>%
  as.matrix()

logfc_mat_all <- logfc_mat_all[, all_contrasts_heatmap, drop = FALSE]
sig_mat_all <- sig_mat_all[, all_contrasts_heatmap, drop = FALSE]
colnames(logfc_mat_all) <- gsub("_", " ", colnames(logfc_mat_all))
colnames(sig_mat_all) <- gsub("_", " ", colnames(sig_mat_all))

if (any(is.na(logfc_mat_all))) {
  stop("NA in logFC matrix for 5-contrast heatmap; check source all_terms CSVs.")
}

# Force exact row order from the CSV exported from rotate heatmap.
target_row_order <- unique(row_annot_tbl$term_label)
missing_rows <- setdiff(target_row_order, rownames(logfc_mat_all))
if (length(missing_rows) > 0L) {
  stop("These term_label rows are not found in matrices: ", paste(missing_rows, collapse = ", "))
}
target_row_order <- target_row_order[target_row_order %in% rownames(logfc_mat_all)]

logfc_mat_all <- logfc_mat_all[target_row_order, , drop = FALSE]
sig_mat_all <- sig_mat_all[target_row_order, , drop = FALSE]
row_annot_tbl <- row_annot_tbl %>%
  filter(term_label %in% target_row_order) %>%
  mutate(term_label = factor(term_label, levels = target_row_order)) %>%
  arrange(term_label) %>%
  mutate(term_label = as.character(term_label))

# Build a rotated dendrogram whose leaf order exactly equals target_row_order.
hc_row_all <- stats::hclust(stats::dist(logfc_mat_all))
d_row <- as.dendrogram(hc_row_all)
# Build cutree blocks first, then rotate blocks by CSV Cluster order (similar to 03 script style).
grp <- stats::cutree(hc_row_all, k = rotate_row_k)
leaf_in_dend_order <- rownames(logfc_mat_all)[hc_row_all$order]
leaf_cluster <- setNames(as.character(row_annot_tbl$cluster), row_annot_tbl$term_label)
cluster_order_csv <- unique(as.character(row_annot_tbl$cluster))
if (nzchar(cluster_last) && cluster_last %in% cluster_order_csv) {
  cluster_order_csv <- c(cluster_order_csv[cluster_order_csv != cluster_last], cluster_last)
}
grp_ids <- sort(unique(as.integer(grp)))
grp_meta <- lapply(grp_ids, function(gid) {
  leaves <- leaf_in_dend_order[grp[leaf_in_dend_order] == gid]
  cls <- leaf_cluster[leaves]
  cls <- cls[!is.na(cls)]
  if (length(cls) < 1L) {
    dom <- NA_character_
  } else {
    tab <- sort(table(cls), decreasing = TRUE)
    dom <- names(tab)[1L]
  }
  first_pos <- min(match(leaves, leaf_in_dend_order))
  dom_rank <- match(dom, cluster_order_csv)
  if (is.na(dom_rank)) dom_rank <- length(cluster_order_csv) + 1L
  tibble::tibble(
    grp_id = gid,
    dominant_cluster = dom,
    dominant_rank = dom_rank,
    first_pos = first_pos
  )
}) %>% dplyr::bind_rows() %>%
  dplyr::arrange(first_pos)
# Keep original block order, only move the requested cluster to bottom.
ordered_grp_ids <- grp_meta$grp_id
if (nzchar(cluster_last) && cluster_last %in% grp_meta$dominant_cluster) {
  move_ids <- grp_meta$grp_id[grp_meta$dominant_cluster == cluster_last]
  stay_ids <- grp_meta$grp_id[grp_meta$dominant_cluster != cluster_last]
  ordered_grp_ids <- c(stay_ids, move_ids)
}
new_leaf_order <- unlist(lapply(ordered_grp_ids, function(gid) {
  leaf_in_dend_order[grp[leaf_in_dend_order] == gid]
}), use.names = FALSE)
d_row_rot <- dendextend::rotate(d_row, order = new_leaf_order)

# cluster: from annotation CSV column Cluster (defines Function and heatmap body gaps; not hclust cutree).

cluster_function_tbl <- row_annot_tbl %>%
  group_by(cluster) %>%
  arrange(is.na(IC), IC, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    cluster,
    Function = term_label,
    Function_term = term,
    Function_GO_id = GO_id,
    Function_IC = IC
  )

if ("Function" %in% colnames(row_annot_tbl)) {
  row_annot_tbl <- row_annot_tbl %>%
    left_join(
      cluster_function_tbl %>%
        rename(
          Function_new = Function,
          Function_term_new = Function_term,
          Function_GO_id_new = Function_GO_id,
          Function_IC_new = Function_IC
        ),
      by = "cluster"
    ) %>%
    mutate(
      Function = dplyr::coalesce(Function, Function_new),
      Function_term = dplyr::coalesce(Function_term, Function_term_new),
      Function_GO_id = dplyr::coalesce(Function_GO_id, Function_GO_id_new),
      Function_IC = dplyr::coalesce(Function_IC, Function_IC_new)
    ) %>%
    dplyr::select(-Function_new, -Function_term_new, -Function_GO_id_new, -Function_IC_new)
} else {
  row_annot_tbl <- row_annot_tbl %>% left_join(cluster_function_tbl, by = "cluster")
}

# GOSemSim BP IC; per row max IC over {term GO ∪ BP ancestors} (export column).
# MaxAncIC PDF: same rows as the %-threshold common-ancestor table: per cluster, take
# the ancestor row with max ancestor_IC (then ancestor_term); ties broken by ancestor_id.
go_sim_ic <- GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = "BP", computeIC = TRUE)
bp_anc_map <- as.list(GO.db::GOBPANCESTOR)
ic_slot <- go_sim_ic@IC
row_annot_tbl$max_ancestor_IC <- vapply(
  seq_len(nrow(row_annot_tbl)),
  function(i) {
    gid <- as.character(row_annot_tbl$GO_id[i])
    if (is.na(gid) || !nzchar(gid) || !grepl("^GO:", gid)) {
      return(NA_real_)
    }
    anc <- unique(c(gid, unname(bp_anc_map[[gid]] %||% character())))
    vals <- as.numeric(unname(ic_slot[anc]))
    vals <- vals[is.finite(vals)]
    if (length(vals) < 1L) {
      return(NA_real_)
    }
    max(vals, na.rm = TRUE)
  },
  numeric(1)
)
drop_broad_bp_ancestor_ids <- c("all", "GO:0008150")
go_bp_term_lookup <- function(go_ids) {
  go_ids <- unique(stats::na.omit(as.character(go_ids)))
  if (length(go_ids) == 0L) {
    return(tibble::tibble(GOID = character(), TERM = character()))
  }
  AnnotationDbi::select(
    GO.db::GO.db,
    keys = go_ids,
    columns = "TERM",
    keytype = "GOID"
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::distinct(GOID, .keep_all = TRUE)
}
ancestor_filter_rows <- list()
for (cl in sort(unique(row_annot_tbl$cluster))) {
  cl_tbl <- row_annot_tbl %>% dplyr::filter(cluster == cl, !is.na(GO_id), grepl("^GO:", GO_id))
  n_cluster_go <- nrow(cl_tbl)
  if (n_cluster_go < 1L) next
  min_n <- ceiling(common_ancestor_frac * n_cluster_go)
  all_pairs <- list()
  for (i in seq_len(nrow(cl_tbl))) {
    gid <- as.character(cl_tbl$GO_id[i])
    anc <- unique(c(gid, unname(bp_anc_map[[gid]] %||% character())))
    anc <- anc[!is.na(anc) & nzchar(anc)]
    for (a in anc) {
      all_pairs[[length(all_pairs) + 1L]] <- tibble::tibble(
        ancestor_id = a,
        term_GO_id = gid,
        term_label = as.character(cl_tbl$term_label[i])
      )
    }
  }
  pairs_tbl <- dplyr::bind_rows(all_pairs)
  if (nrow(pairs_tbl) < 1L) next
  go_term_map <- go_bp_term_lookup(
    c(unique(pairs_tbl$term_GO_id), unique(pairs_tbl$ancestor_id))
  )
  id2term <- setNames(go_term_map$TERM, go_term_map$GOID)
  out_cl <- pairs_tbl %>%
    dplyr::filter(!as.character(ancestor_id) %in% drop_broad_bp_ancestor_ids) %>%
    dplyr::count(ancestor_id, name = "n_terms_with_ancestor", sort = TRUE) %>%
    dplyr::filter(n_terms_with_ancestor >= min_n) %>%
    dplyr::mutate(
      cluster = cl,
      n_terms_cluster = n_cluster_go,
      min_terms_required = min_n,
      frac = n_terms_with_ancestor / n_terms_cluster,
      ancestor_term = unname(id2term[as.character(ancestor_id)]),
      ancestor_IC = as.numeric(unname(go_sim_ic@IC[as.character(ancestor_id)]))
    ) %>%
    dplyr::left_join(
      pairs_tbl %>%
        dplyr::group_by(ancestor_id) %>%
        dplyr::summarise(
          term_GO_ids = paste(sort(unique(term_GO_id)), collapse = "; "),
          term_labels = paste(sort(unique(term_label)), collapse = "; "),
          .groups = "drop"
        ),
      by = "ancestor_id"
    ) %>%
    dplyr::arrange(dplyr::desc(frac), ancestor_id) %>%
    dplyr::select(
      cluster,
      n_terms_cluster,
      min_terms_required,
      ancestor_id,
      ancestor_term,
      ancestor_IC,
      n_terms_with_ancestor,
      frac,
      term_GO_ids,
      term_labels
    )
  ancestor_filter_rows[[length(ancestor_filter_rows) + 1L]] <- out_cl
}
cluster_bp_ancestor_by_pct <- dplyr::bind_rows(ancestor_filter_rows)
if (nrow(cluster_bp_ancestor_by_pct) >= 1L) {
  rep_from_anc_table <- cluster_bp_ancestor_by_pct %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::arrange(dplyr::desc(.data$ancestor_IC), .data$ancestor_id) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      cluster = .data$cluster,
      MaxAncIC_ancestor_id = .data$ancestor_id,
      MaxAncIC_ancestor_term = .data$ancestor_term,
      Cluster_winner_ancestor_IC = .data$ancestor_IC
    )
} else {
  rep_from_anc_table <- tibble::tibble(
    cluster = character(0L),
    MaxAncIC_ancestor_id = character(0L),
    MaxAncIC_ancestor_term = character(0L),
    Cluster_winner_ancestor_IC = double(0L)
  )
}

# Rotate order for *_rotate heatmaps: enforce MaxAncIC block order left->right.
desired_maxanc_order <- c(
  "regulation of immune system process",
  "biological regulation",
  "metabolic process",
  "protein maturation",
  "extrinsic apoptotic signaling pathway",
  "leukocyte chemotaxis"
)
if (nrow(rep_from_anc_table) >= 1L) {
  cl_term <- rep_from_anc_table %>%
    dplyr::mutate(maxanc_term_lc = tolower(as.character(.data$MaxAncIC_ancestor_term)))
  clusters_in_target_order <- unlist(lapply(desired_maxanc_order, function(tt) {
    cl_term$cluster[cl_term$maxanc_term_lc == tt]
  }), use.names = FALSE)
  clusters_in_target_order <- unique(clusters_in_target_order)
  cluster_order_rotate <- c(
    clusters_in_target_order,
    cluster_order_csv[!cluster_order_csv %in% clusters_in_target_order]
  )
  ordered_grp_ids <- unlist(lapply(cluster_order_rotate, function(cl) {
    grp_meta$grp_id[grp_meta$dominant_cluster == cl]
  }), use.names = FALSE)
  ordered_grp_ids <- c(
    ordered_grp_ids,
    grp_meta$grp_id[!grp_meta$grp_id %in% ordered_grp_ids]
  )
  ordered_grp_ids <- unique(ordered_grp_ids)
  new_leaf_order <- unlist(lapply(ordered_grp_ids, function(gid) {
    leaf_in_dend_order[grp[leaf_in_dend_order] == gid]
  }), use.names = FALSE)
  d_row_rot <- dendextend::rotate(d_row, order = new_leaf_order)
}
row_annot_tbl <- row_annot_tbl %>% dplyr::left_join(rep_from_anc_table, by = "cluster")
# Legend colors by displayed string (unresolved GO -> placeholder).
row_annot_tbl$MaxAncIC_ancestor_term_disp <- as.character(row_annot_tbl$MaxAncIC_ancestor_term)
miss <- is.na(row_annot_tbl$MaxAncIC_ancestor_term_disp) | !nzchar(row_annot_tbl$MaxAncIC_ancestor_term_disp)
row_annot_tbl$MaxAncIC_ancestor_term_disp[miss] <- "\u2014"
maxanc_levels <- unique(row_annot_tbl$MaxAncIC_ancestor_term_disp)
if (length(maxanc_levels) > length(my_palette)) {
  stop("my_palette does not have enough colors for MaxAncIC ancestor TERM labels.")
}
maxanc_cols <- setNames(my_palette[seq_along(maxanc_levels)], maxanc_levels)
row_annot_tbl$MaxAncIC_ancestor_term_disp <- factor(
  row_annot_tbl$MaxAncIC_ancestor_term_disp,
  levels = maxanc_levels
)
ha_left_maxanc <- ComplexHeatmap::rowAnnotation(
  MaxAncIC = row_annot_tbl$MaxAncIC_ancestor_term_disp,
  col = list(MaxAncIC = maxanc_cols),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  width = unit(6, "mm")
)

fc_range_all <- range(logfc_mat_all, na.rm = TRUE)
fc_max_all <- max(abs(fc_range_all))
col_logfc_all <- circlize::colorRamp2(
  c(-fc_max_all, 0, fc_max_all),
  c("navy", "white", "firebrick3")
)

fun_levels <- unique(row_annot_tbl$Function)
if (length(fun_levels) > length(my_palette)) {
  stop("my_palette does not have enough colors for Function annotations.")
}
fun_cols <- setNames(my_palette[seq_along(fun_levels)], fun_levels)
fun_by_term <- setNames(as.character(row_annot_tbl$Function), row_annot_tbl$term_label)
maxanc_by_term <- setNames(as.character(row_annot_tbl$MaxAncIC_ancestor_term_disp), row_annot_tbl$term_label)
maxanc_cols <- setNames(my_palette[seq_along(maxanc_levels)], maxanc_levels)

# Horizontal layout: contrasts on rows, terms on columns
logfc_mat_plot <- t(logfc_mat_all)
sig_mat_plot <- t(sig_mat_all)

ha_top_fun <- ComplexHeatmap::HeatmapAnnotation(
  Function = factor(fun_by_term[colnames(logfc_mat_plot)], levels = fun_levels),
  col = list(Function = fun_cols),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)
ha_top_maxanc <- ComplexHeatmap::HeatmapAnnotation(
  MaxAncIC = factor(maxanc_by_term[colnames(logfc_mat_plot)], levels = maxanc_levels),
  col = list(MaxAncIC = maxanc_cols),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

ht_all_fun <- ComplexHeatmap::Heatmap(
  logfc_mat_plot,
  name = "logFC",
  col = col_logfc_all,
  top_annotation = ha_top_fun,
  width = unit(ncol(logfc_mat_plot) * 6, "mm"),
  height = unit(nrow(logfc_mat_plot) * 8, "mm"),
  cluster_rows = FALSE,
  cluster_columns = d_row,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 90,
  cell_fun = function(j, i, x, y, width, height, fill) {
    lab <- sig_mat_plot[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      grid.text(lab, x, y, gp = gpar(fontsize = 8, fontface = "bold", col = "black"))
    }
  },
  heatmap_legend_param = list(
    title = "logFC",
    direction = "vertical"
  )
)

ht_all_rot_fun <- ComplexHeatmap::Heatmap(
  logfc_mat_plot,
  name = "logFC",
  col = col_logfc_all,
  top_annotation = ha_top_fun,
  width = unit(ncol(logfc_mat_plot) * 6, "mm"),
  height = unit(nrow(logfc_mat_plot) * 8, "mm"),
  cluster_rows = FALSE,
  cluster_columns = d_row_rot,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 90,
  cell_fun = function(j, i, x, y, width, height, fill) {
    lab <- sig_mat_plot[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      grid.text(lab, x, y, gp = gpar(fontsize = 8, fontface = "bold", col = "black"))
    }
  },
  heatmap_legend_param = list(
    title = "logFC",
    direction = "vertical"
  )
)

pdf(
  file.path(
    output_dir,
    "GO_BP_ssGSEA_top10_5comparisons_logFC_with_stars.pdf"
  ),
  width = 14,
  height = max(6, nrow(logfc_mat_all) * 0.26)
)
ComplexHeatmap::draw(
  ht_all_fun,
  annotation_legend_side = "right",
  heatmap_legend_side = "left",
  padding = unit(c(2, 2, 2, 2), "mm")
)
# White gaps between cluster blocks in displayed column order.
cluster_by_term <- setNames(as.character(row_annot_tbl$cluster), row_annot_tbl$term_label)
order_nonrot <- rownames(logfc_mat_all)[stats::order.dendrogram(d_row)]
cl_vec_nonrot <- unname(cluster_by_term[order_nonrot])
cluster_boundary_idx_nonrot <- which(cl_vec_nonrot[-1L] != cl_vec_nonrot[-length(cl_vec_nonrot)])
if (length(cluster_boundary_idx_nonrot) > 0L) {
  ComplexHeatmap::decorate_heatmap_body("logFC", {
    n <- ncol(logfc_mat_plot)
    for (b in cluster_boundary_idx_nonrot) {
      x <- b / n
      grid.lines(
        x = unit(c(x, x), "npc"),
        y = unit(c(0, 1), "npc"),
        gp = gpar(col = "white", lwd = 5)
      )
    }
  })
}
dev.off()

# Additional cluster-rotate version (same rotated row dendrogram, explicit filename).
pdf(
  file.path(
    output_dir,
    "GO_BP_ssGSEA_top10_5comparisons_logFC_with_stars_rowCluster_rotate.pdf"
  ),
  width = 14,
  height = max(6, nrow(logfc_mat_all) * 0.26)
)
ComplexHeatmap::draw(
  ht_all_rot_fun,
  annotation_legend_side = "right",
  heatmap_legend_side = "left",
  padding = unit(c(2, 2, 2, 2), "mm")
)
order_rot <- rownames(logfc_mat_all)[stats::order.dendrogram(d_row_rot)]
cl_vec_rot <- unname(cluster_by_term[order_rot])
cluster_boundary_idx_rot <- which(cl_vec_rot[-1L] != cl_vec_rot[-length(cl_vec_rot)])
if (length(cluster_boundary_idx_rot) > 0L) {
  ComplexHeatmap::decorate_heatmap_body("logFC", {
    n <- ncol(logfc_mat_plot)
    for (b in cluster_boundary_idx_rot) {
      x <- b / n
      grid.lines(
        x = unit(c(x, x), "npc"),
        y = unit(c(0, 1), "npc"),
        gp = gpar(col = "white", lwd = 5)
      )
    }
  })
}
dev.off()

# Same heatmap; left = ancestor_term with max ancestor_IC within that cluster in the %%-threshold
# common-ancestor table (same as BP_ancestors_at_least_*pct CSV).
pdf(
  file.path(
    output_dir,
    "GO_BP_ssGSEA_top10_5comparisons_logFC_with_stars_annotMaxAncICterm.pdf"
  ),
  width = 14,
  height = max(6, nrow(logfc_mat_all) * 0.26)
)
ht_all_fun_maxanc <- ComplexHeatmap::Heatmap(
  logfc_mat_plot,
  name = "logFC",
  col = col_logfc_all,
  top_annotation = ha_top_maxanc,
  width = unit(ncol(logfc_mat_plot) * 6, "mm"),
  height = unit(nrow(logfc_mat_plot) * 8, "mm"),
  cluster_rows = FALSE,
  cluster_columns = d_row,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 90,
  cell_fun = function(j, i, x, y, width, height, fill) {
    lab <- sig_mat_plot[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      grid.text(lab, x, y, gp = gpar(fontsize = 8, fontface = "bold", col = "black"))
    }
  },
  heatmap_legend_param = list(
    title = "logFC",
    direction = "vertical"
  )
)
ComplexHeatmap::draw(
  ht_all_fun_maxanc,
  annotation_legend_side = "right",
  heatmap_legend_side = "left",
  padding = unit(c(2, 2, 2, 2), "mm")
)
if (length(cluster_boundary_idx_nonrot) > 0L) {
  ComplexHeatmap::decorate_heatmap_body("logFC", {
    n <- ncol(logfc_mat_plot)
    for (b in cluster_boundary_idx_nonrot) {
      x <- b / n
      grid.lines(
        x = unit(c(x, x), "npc"),
        y = unit(c(0, 1), "npc"),
        gp = gpar(col = "white", lwd = 5)
      )
    }
  })
}
dev.off()

# Additional cluster-rotate version for MaxAncIC annotation.
pdf(
  file.path(
    output_dir,
    "GO_BP_ssGSEA_top10_5comparisons_logFC_with_stars_annotMaxAncICterm_rowCluster_rotate.pdf"
  ),
  width = 14,
  height = max(6, nrow(logfc_mat_all) * 0.26)
)
ht_all_rot_fun_maxanc <- ComplexHeatmap::Heatmap(
  logfc_mat_plot,
  name = "logFC",
  col = col_logfc_all,
  top_annotation = ha_top_maxanc,
  width = unit(ncol(logfc_mat_plot) * 6, "mm"),
  height = unit(nrow(logfc_mat_plot) * 8, "mm"),
  cluster_rows = FALSE,
  cluster_columns = d_row_rot,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 90,
  cell_fun = function(j, i, x, y, width, height, fill) {
    lab <- sig_mat_plot[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      grid.text(lab, x, y, gp = gpar(fontsize = 8, fontface = "bold", col = "black"))
    }
  },
  heatmap_legend_param = list(
    title = "logFC",
    direction = "vertical"
  )
)
ComplexHeatmap::draw(
  ht_all_rot_fun_maxanc,
  annotation_legend_side = "right",
  heatmap_legend_side = "left",
  padding = unit(c(2, 2, 2, 2), "mm")
)
if (length(cluster_boundary_idx_rot) > 0L) {
  ComplexHeatmap::decorate_heatmap_body("logFC", {
    n <- ncol(logfc_mat_plot)
    for (b in cluster_boundary_idx_rot) {
      x <- b / n
      grid.lines(
        x = unit(c(x, x), "npc"),
        y = unit(c(0, 1), "npc"),
        gp = gpar(col = "white", lwd = 5)
      )
    }
  })
}
dev.off()

# Export: same columns as rotate_cluster_IC_per_term + Function, Function_GO_id (row order = heatmap).
# cluster_bp_ancestor_by_pct: already built (before heatmaps) for MaxAncIC + this CSV.
if ("Cluster" %in% names(row_annot_tbl)) {
  row_annot_export <- row_annot_tbl %>%
    dplyr::select(
      "term_label",
      "term",
      "GO_id",
      "IC",
      "Cluster",
      "Function",
      "Function_GO_id"
    )
} else {
  row_annot_export <- row_annot_tbl %>%
    dplyr::transmute(
      .data$term_label,
      .data$term,
      .data$GO_id,
      .data$IC,
      Cluster = as.character(.data$cluster),
      .data$Function,
      .data$Function_GO_id
    )
}
readr::write_csv(
  row_annot_export,
  file.path(output_dir, "GO_BP_ssGSEA_top10_5comparisons_rotate_withFunction_rowOrdered.csv")
)
readr::write_csv(
  cluster_bp_ancestor_by_pct,
  file.path(
    output_dir,
    sprintf(
      "GO_BP_ssGSEA_top10_5comparisons_rotate_cluster_BP_ancestors_at_least_%spct_terms.csv",
      pct_fname_tag
    )
  )
)

message(
  "Common-ancestor table: >= ",
  common_ancestor_pct,
  "% of terms per cluster (min count = ceiling(fraction * n_terms))."
)
message("Saved outputs to: ", output_dir)
