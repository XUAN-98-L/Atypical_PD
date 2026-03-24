gsva_results = readRDS("Results/00_aPS_Ctrl_PD_paired_comparison/BP_ssGSEA_results.rds")

library(limma)
library(DEP)
library(tidyverse)
# Define the Z-score normalization function
zscore_normalize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Apply the Z-score normalization to each row (gene set)
gsva_results <- t(apply(gsva_results, 1, zscore_normalize))

metadata = readxl::read_excel("Data/CSF_APS_individual_patient_results.xlsx",sheet = 2)

group_info = metadata %>% dplyr::select("ID","Subtype")
# replace "Non-PD-Synucleinopathy" into "NPS"
group_info$Subtype[group_info$Subtype == "Non-PD-Synucleinopathy"] = "NPS"
# exclude Ctr from modeling
group_info <- group_info %>% dplyr::filter(Subtype != "Ctr")
common_ids <- intersect(colnames(gsva_results), group_info$ID)
if (length(common_ids) == 0) {
  stop("No overlapping sample IDs between gsva_results and metadata after excluding Ctr.")
}
gsva_results <- gsva_results[, common_ids, drop = FALSE]
group_info <- group_info %>% dplyr::filter(ID %in% common_ids)
group_info <- group_info[match(colnames(gsva_results), group_info$ID), ]

# Apply limma (NPS, PD, Tauopathy only)
group_info$Subtype <- droplevels(as.factor(group_info$Subtype))
limma_design <- model.matrix(~ 0 + group_info$Subtype)
colnames(limma_design) <- levels(group_info$Subtype)
# Fit the linear model
fit <- lmFit(gsva_results, limma_design)
# Contrasts without Ctr: disease vs disease only
contrast_matrix <- makeContrasts(NPS-PD, Tauopathy-PD, Tauopathy-NPS, levels = limma_design)

# Apply contrasts to the fit
fit2 <- contrasts.fit(fit, contrast_matrix)
# Compute empirical Bayes statistics
fit2 <- eBayes(fit2)

# Extract top pathways per contrast (p < 0.1 for heatmap); no Ctr contrasts
NPS_vs_PD <- topTable(fit2, coef = "NPS - PD", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE, p.value = 0.1) %>% rownames_to_column("term")
Tauopathy_vs_PD <- topTable(fit2, coef = "Tauopathy - PD", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE, p.value = 0.1) %>% rownames_to_column("term")
Tauopathy_vs_NPS <- topTable(fit2, coef = "Tauopathy - NPS", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE, p.value = 0.1) %>% rownames_to_column("term")

# Full tables (no p-value cutoff) for export
NPS_vs_PD_all <- topTable(fit2, coef = "NPS - PD", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE) %>% rownames_to_column("term")
Tauopathy_vs_PD_all <- topTable(fit2, coef = "Tauopathy - PD", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE) %>% rownames_to_column("term")
Tauopathy_vs_NPS_all <- topTable(fit2, coef = "Tauopathy - NPS", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE) %>% rownames_to_column("term")


# Load dplyr package
library(dplyr)
library(ComplexHeatmap)
library(circlize)

############ Prepare data for heatmap
# Split into top and bottom halves based on t values
topn = 10
NPS_vs_PD_up = NPS_vs_PD %>% dplyr::arrange(dplyr::desc(.data[["t"]])) %>% dplyr::slice_head(n = topn / 2) %>% as.data.frame()
NPS_vs_PD_down = NPS_vs_PD %>% dplyr::arrange(.data[["t"]]) %>% dplyr::slice_head(n = topn / 2) %>% as.data.frame()
NPS_vs_PD_forplot = rbind(NPS_vs_PD_up,NPS_vs_PD_down)%>% dplyr::distinct()

Tauopathy_vs_PD_up = Tauopathy_vs_PD %>% dplyr::arrange(dplyr::desc(.data[["t"]])) %>% dplyr::slice_head(n = topn / 2) %>% as.data.frame()
Tauopathy_vs_PD_down = Tauopathy_vs_PD %>% dplyr::arrange(.data[["t"]]) %>% dplyr::slice_head(n = topn / 2) %>% as.data.frame()
Tauopathy_vs_PD_forplot = rbind(Tauopathy_vs_PD_up,Tauopathy_vs_PD_down)%>% dplyr::distinct()

# Tauopathy_vs_NPS_up = Tauopathy_vs_NPS %>% dplyr::arrange(dplyr::desc(.data[["t"]])) %>% dplyr::slice_head(n = topn / 2) %>% as.data.frame()
# Tauopathy_vs_NPS_down = Tauopathy_vs_NPS %>% dplyr::arrange(.data[["t"]]) %>% dplyr::slice_head(n = topn / 2) %>% as.data.frame()
# Tauopathy_vs_NPS_forplot = rbind(Tauopathy_vs_NPS_up,Tauopathy_vs_NPS_down)%>% dplyr::distinct()

topn_terms = rbind(NPS_vs_PD_forplot %>% mutate(contrast = "NPS_vs_PD"),
                   Tauopathy_vs_PD_forplot %>% mutate(contrast = "Tauopathy_vs_PD")
                   #Tauopathy_vs_NPS_forplot %>% mutate(contrast = "Tauopathy_vs_NPS")
                   )

output_dir = "Results/01_Top10_6comparisons_without_CTRL"
if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}
sub_name = "GO_BP"
write.csv(topn_terms, quote = F, row.names = F, file = file.path(output_dir, paste0(sub_name, "_ssGSEA_top10_terms_for_heatmap.csv")))

# Export full tables (all terms, no p-value cutoff) for each comparison
write.csv(NPS_vs_PD_all, quote = F, row.names = F, file = file.path(output_dir, paste0(sub_name, "_ssGSEA_NPS_vs_PD_all_terms.csv")))
write.csv(Tauopathy_vs_PD_all, quote = F, row.names = F, file = file.path(output_dir, paste0(sub_name, "_ssGSEA_Tauopathy_vs_PD_all_terms.csv")))
write.csv(Tauopathy_vs_NPS_all, quote = F, row.names = F, file = file.path(output_dir, paste0(sub_name, "_ssGSEA_Tauopathy_vs_NPS_all_terms.csv")))

#rownames(gsva_results) %in% topn_terms$term

# Sample-level heatmap (all samples are non-Ctr after filtering above)
annolabel <- data.frame(Group = group_info[, "Subtype"])
rownames(annolabel) <- group_info$ID
ordered_columns <- rownames(annolabel)[order(annolabel$Subtype)]
group_colors <- list(Subtype = c(
  "NPS" = "#8DA0CB",
  "PD" = "#A6D854",
  "Tauopathy" = "#E78AC3"
))

plot_df <- gsva_results[
  unique(c(NPS_vs_PD_forplot$term, Tauopathy_vs_PD_forplot$term)), #, Tauopathy_vs_NPS_forplot$term)),
  ordered_columns
]


# Custom function to abbreviate long terms by adding "..." at the end
abbreviate_terms_end <- function(term, max_length = 60) {
  if (nchar(term) > max_length) {
    paste0(substr(term, 1, max_length - 3), "...")
  } else {
    term
  }
}

# Extract the row names from metadata according to the order of colnames in plot_df
rownames(plot_df) <- sapply(rownames(plot_df), abbreviate_terms_end)

#sample level ssGSEA heatmap
pathway_heatmap <- pheatmap::pheatmap(plot_df,
                                      annotation_col = annolabel,
                                      annotation_colors = group_colors,  
                                      fontsize_row = 6,fontsize_col = 6, cluster_cols = FALSE, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), border_color = NA)

pdf(file.path(output_dir, paste0(sub_name,"_ssGSEA_sample_level.pdf")),width = 12)
print(pathway_heatmap)

### Means heatmap (no Ctr subtype column)
sam_long <- as.data.frame(t(gsva_results))
sam_long$ID <- rownames(sam_long)
sam_long <- dplyr::left_join(sam_long, group_info, by = "ID")
tmp_means <- sam_long %>%
  dplyr::select(-ID) %>%
  dplyr::group_by(Subtype) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), mean, na.rm = TRUE)) %>%
  as.data.frame()
rownames(tmp_means) <- tmp_means$Subtype
tmp_means <- tmp_means[, colnames(tmp_means) != "Subtype", drop = FALSE]
Subtype_means <- t(as.matrix(tmp_means))

write_Subtype_means = Subtype_means %>% as.data.frame() %>%rownames_to_column("terms")
write.csv(write_Subtype_means,quote = F,row.names = F, file = file.path(output_dir, paste0(sub_name,"_ssGSEA_Subtype_means.csv")))

Subtype_means = Subtype_means[unique(c(NPS_vs_PD_forplot$term, Tauopathy_vs_PD_forplot$term)), ] #, Tauopathy_vs_NPS_forplot$term)),

# Custom function to abbreviate long terms by adding "..." at the end
rownames(Subtype_means) <- sapply(rownames(Subtype_means), abbreviate_terms_end)
# remove the "GOBP_" prefix
rownames(Subtype_means) <- gsub("GOBP_", "", rownames(Subtype_means))
# remove "_" and replace with " "
rownames(Subtype_means) <- gsub("_", " ", rownames(Subtype_means))
# change into lower case
rownames(Subtype_means) <- tolower(rownames(Subtype_means))
# change into only the first letter capitalized
rownames(Subtype_means) <- tools::toTitleCase(rownames(Subtype_means))

# Column annotation for Subtype means heatmap
subtype_colors <- c("NPS" = "#4975FF", "PD" = "#E9967B", "Tauopathy" = "#BF3EFF")
ha_mean <- HeatmapAnnotation(Subtype = colnames(Subtype_means), col = list(Subtype = subtype_colors), show_annotation_name = FALSE)

# Continuous color scale for z-scores (same as navy/white/firebrick3)
z_range <- range(Subtype_means, na.rm = TRUE)
col_means <- circlize::colorRamp2(c(z_range[1], 0, z_range[2]), c("navy", "white", "firebrick3"))

pathway_heatmap_mean <- ComplexHeatmap::Heatmap(Subtype_means,
  name = "Mean Z-scored ssGSEA scores",
  col = col_means,
  top_annotation = ha_mean,
  cluster_columns = FALSE,
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 8),
  heatmap_legend_param = list(
    title = "Mean Z-scored ssGSEA scores",
    title_position = "leftcenter-rot",
    legend_direction = "vertical"
  )
)

pdf(file.path(output_dir, paste0(sub_name, "_ssGSEA_Subtype_means.pdf")), width = 7)
draw(pathway_heatmap_mean, heatmap_legend_side = "left")
dev.off()
