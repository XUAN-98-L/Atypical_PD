library(tidyverse)
library(msigdbr)
library(biomaRt)
library(limma)
library(writexl)
library(GSVA)

# Loading data
expression_data = readxl::read_excel("Data/CSF_APS_individual_patient_results.xlsx",sheet = 1)
metadata = readxl::read_excel("Data/CSF_APS_individual_patient_results.xlsx",sheet = 2)
output_dir = "Results/00_aPS_Ctrl_PD_paired_comparison"
if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

group_info = metadata %>% dplyr::select("ID","Condition")

##########  protein name cleaning
protein_expression <- expression_data %>%
  mutate(id_first = gsub(";.*", "", `ID`)) %>%
  tidyr::separate(id_first, into = c("SwissProt_Tremble", "UniProtAccession", "UniProtName"), sep = "\\|") %>%
  dplyr::filter(SwissProt_Tremble == "sp", UniProtName != "Biognosys_iRT")

protein_expression = protein_expression %>% dplyr::select(-c("ID", "SwissProt_Tremble","UniProtName")) %>% column_to_rownames(var = "UniProtAccession")

########## protein name to gene symbol
# Download https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz
# HUMAN_9606_idmapping_selected.tab.gz	2024-11-27 14:00	119M	
database = "Data/HUMAN_9606_idmapping_selected.tab.gz"
message("Reading protein to gene mapping...")
prot_to = readr::read_tsv(database, col_names = FALSE)

prot_to_gene <- prot_to[,c("X1","X2", "X19")]
colnames(prot_to_gene) <- c("UniProtAccession","UniProtName","gene_id")

# Find intersecting UniProt accessions between 'protein_expression' and 'prot_to_gene'
inter = intersect(rownames(protein_expression), prot_to_gene$UniProtAccession)
# Filter 'prot_to_gene' to include only the matching UniProt accessions found in 'protein_expression'
match_ids <- prot_to_gene[which(prot_to_gene$UniProtAccession %in% inter),]
# Sort 'match_ids' by 'UniProtAccession' to ensure it matches the order in 'protein_expression'
match_ids = match_ids[match(rownames(protein_expression),match_ids$UniProtAccession),]
# Remove the "_HUMAN" suffix from 'UniProtName' to obtain a clean 'gene_name'
match_ids$gene_name = gsub("_HUMAN","",match_ids$UniProtName)
# define gene ID use for search
gene_ids <- match_ids$UniProtAccession
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_info = getBM(
  attributes=c('ensembl_gene_id','uniprotswissprot','hgnc_symbol'), 
  mart = mart) 

# rename specific columns
colnames(gene_info)[colnames(gene_info) == "hgnc_symbol"] <- "SYMBOL"
colnames(gene_info)[colnames(gene_info) == "uniprotswissprot"] <- "UniProtAccession"

# keep only unique uniprot accession numbers
gene_info = gene_info %>% distinct(gene_info$UniProtAccession, .keep_all = TRUE)

df = left_join(match_ids, gene_info, by = "UniProtAccession" )

# if df$SYMBOL is "", replace it into NA
df$SYMBOL[df$SYMBOL == ""] = NA
# if df$SYMBOL is not NA, use SYMBOL instead of gene_name
df$gene_name[!is.na(df$SYMBOL)] = df$SYMBOL[!is.na(df$SYMBOL)]

#rename SYMBOL in df into better_gene_name
uniprot_to_genename <- df %>%
  dplyr::rename(better_gene_name = SYMBOL) %>%
  dplyr::select(ensembl_gene_id,gene_name,UniProtAccession,better_gene_name,UniProtName)

colnames(uniprot_to_genename)[colnames(uniprot_to_genename) == "ensembl_gene_id"] <- "gene_id"

# make rows in uniprot_to_genename the same order as in protein_expression matrix
uniprot_to_genename = uniprot_to_genename[match(rownames(protein_expression),uniprot_to_genename$UniProtAccession),]

protein_expression$gene_name = uniprot_to_genename$gene_name
rownames(protein_expression) = protein_expression$gene_name
protein_expression = protein_expression %>% dplyr::select(-c("gene_name"))

protein_expression = protein_expression[,group_info$ID]

###prepare geneset for ssGSEA
category = "C5"
message(paste("msigdbr category used in GSVA is:",category))

for (subset in c("GO:BP","GO:CC","GO:MF")){
  sub_name = str_split(subset,":")[[1]][2]
  GO_gene_sets <- msigdbr::msigdbr(
    species = "Homo sapiens", # Can change this to what species you need
    category = category # Only GO gene sets
  ) 
  GO_gene_sets = GO_gene_sets[which(GO_gene_sets$gs_subcat ==subset),]
  
  # Background removal
  GO_gene_sets <- GO_gene_sets[GO_gene_sets$gene_symbol %in% rownames(protein_expression), ]
  
  GO_list <- split(
    GO_gene_sets$gene_symbol, # The genes we want split into pathways
    GO_gene_sets$gs_name # The pathways made as the higher levels of the list
  )
  
  filtered_mapped_matrix = protein_expression%>%
    mutate_all(as.numeric)%>% as.matrix()
  
  # Extract gene names from filtered_mapped_matrix
  genes_of_interest <- rownames(filtered_mapped_matrix)
  
  # Function to calculate overlap with genes of interest
  overlap_with_genes <- function(geneset) {
    length(intersect(geneset, genes_of_interest))
  }
  
  # Filter GO_list, keeping only sets with at least 10 overlapping genes
  filtered_GO_list <- GO_list[sapply(GO_list, overlap_with_genes) >= 10]
  
  gsva_results <- gsva(
    filtered_mapped_matrix,
    filtered_GO_list,
    method = "ssgsea",
    # Appropriate for our vst transformed data
    kcdf = "Gaussian",
    # Minimum gene set size
    min.sz = 5,
    # Maximum gene set size
    max.sz = 500,
    # Compute Gaussian-distributed scores
    mx.diff = TRUE,
    # Don't print out the progress bar
    verbose = FALSE
  )
  
  # save GSVA result
  saveRDS(gsva_results,file = file.path(output_dir,paste0(sub_name,"_ssGSEA_results.rds")))
  
  # Define the Z-score normalization function
  zscore_normalize <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }
  
  # Apply the Z-score normalization to each row (gene set)
  gsva_results <- t(apply(gsva_results, 1, zscore_normalize))
  
  ####limma comparison
  # Apply limma for all pairwise comparisons
  group_info$Condition = as.factor(group_info$Condition)
  limma_design <- model.matrix(~ 0 + group_info$Condition)
  colnames(limma_design) <- levels(group_info$Condition)
  # Fit the linear model
  fit <- lmFit(gsva_results, limma_design)
  # Create contrasts for group comparisons
  contrast_matrix <- makeContrasts(aPS-Ctrl, PD-Ctrl,aPS-PD, levels = limma_design)
  # Apply contrasts to the fit
  fit2 <- contrasts.fit(fit, contrast_matrix)
  # Compute empirical Bayes statistics
  fit2 <- eBayes(fit2)
  
  # Extract the top differentially expressed pathways for each contrast
  aPS_vs_Ctrl <- topTable(fit2, coef = "aPS - Ctrl", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE,p.value=0.1)
  PD_vs_Ctrl <- topTable(fit2, coef = "PD - Ctrl", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE,p.value=0.1)
  aPS_vs_PD <- topTable(fit2, coef = "aPS - PD", sort.by = "t", adjust.method = "BH", number = Inf, confint = TRUE,p.value=0.1)
  
  # Define a function to process each comparison dataframe
  process_comparison <- function(comparison_df, comparison_name, GO_list, filtered_mapped_matrix, topn, output_dir) {
    
    # Select top terms from the comparison dataframe and find proteins match to genesets
    comparison_topn_up <- comparison_df %>% dplyr::slice_max(n = (topn+20)/2,order_by = t)  %>% as.data.frame()
    comparison_topn_down <- comparison_df %>% dplyr::slice_min(n = (topn+20)/2,order_by = t)  %>% as.data.frame()
    comparison_topn_df = rbind(comparison_topn_up,comparison_topn_down)%>% dplyr::distinct()
    
    # Find common GO terms between GO_list and gsva_result
    common_GO_terms <- intersect(names(GO_list), rownames(comparison_topn_df))
    
    # Subset GO_list to keep only the terms present in common_GO_terms
    GO_list_subset <- GO_list[names(GO_list) %in% common_GO_terms]
    
    # Extract the rownames of the assay genes
    assay_genes <- rownames(filtered_mapped_matrix)
    
    # Subset genes in GO_list_subset to keep only those present in assay_genes
    GO_list_subset_filtered <- lapply(GO_list_subset, function(genes) {
      intersect(genes, assay_genes)
    })
    
    # Get top terms and genes
    GO_list_top_terms <- GO_list_subset_filtered[names(GO_list_subset_filtered) %in% common_GO_terms]
    
    # Create a data frame with the comparison_df terms and their corresponding genes
    comparison_topn <- comparison_topn_df %>%
      rownames_to_column("terms") %>%
      mutate(genes = map(terms, ~ GO_list_top_terms[[.]])) %>%
      dplyr::select(terms, genes, everything())
    
    # Combine the genes for each term
    comparison_topn$genes <- sapply(comparison_topn$genes, function(x) paste(x, collapse = ","))
    
    # Write the result into an Excel file
    write_xlsx(comparison_topn %>% as.data.frame(),
               path = file.path(output_dir, paste0(sub_name,"_ssGSEA_topn_terms_", comparison_name, "_with_matched_proteins.xlsx")))
    
    return(comparison_topn)
  }
  
  topn = 20
  aPS_vs_Ctrl_topn <- process_comparison(aPS_vs_Ctrl, "aPS_vs_Ctrl", GO_list, filtered_mapped_matrix, topn, output_dir)
  PD_vs_Ctrl_topn <- process_comparison(PD_vs_Ctrl, "PD_vs_Ctrl", GO_list, filtered_mapped_matrix, topn, output_dir)
  aPS_vs_PD_topn <- process_comparison(aPS_vs_PD, "aPS_vs_PD", GO_list, filtered_mapped_matrix, topn, output_dir)
  
  
  ############ Prepare data for heatmap
  # Split aPS_vs_Ctrl into top and bottom halves based on t values
  aPS_vs_Ctrl_up = aPS_vs_Ctrl %>% dplyr::slice_max(n = topn/2,order_by = t)  %>% as.data.frame()
  aPS_vs_Ctrl_down = aPS_vs_Ctrl %>% dplyr::slice_min(n = topn/2,order_by = t)  %>% as.data.frame()
  aPS_vs_Ctrl = rbind(aPS_vs_Ctrl_up,aPS_vs_Ctrl_down)%>% dplyr::distinct()
  
  # Split PD_vs_Ctrl into top and bottom halves based on t values
  PD_vs_Ctrl_up = PD_vs_Ctrl %>% dplyr::slice_max(n = topn/2, order_by = t) %>% as.data.frame()
  PD_vs_Ctrl_down = PD_vs_Ctrl %>% dplyr::slice_min(n = topn/2, order_by = t) %>% as.data.frame()
  #need distinct() to avoid duplicated terms (when there are less enriched terms than topn/2)
  PD_vs_Ctrl = rbind(PD_vs_Ctrl_up, PD_vs_Ctrl_down)%>% dplyr::distinct()
  
  # Split aPS_vs_PD into top and bottom halves based on t values
  aPS_vs_PD_up = aPS_vs_PD %>% dplyr::slice_max(n = topn/2, order_by = t) %>% as.data.frame()
  aPS_vs_PD_down = aPS_vs_PD %>% dplyr::slice_min(n = topn/2, order_by = t) %>% as.data.frame()
  aPS_vs_PD = rbind(aPS_vs_PD_up, aPS_vs_PD_down)%>% dplyr::distinct()
  
  topn_terms = rbind(aPS_vs_Ctrl %>% mutate(contrast = "aPS_vs_Ctrl"),
                     PD_vs_Ctrl %>% mutate(contrast = "PD_vs_Ctrl"),
                     aPS_vs_PD %>% mutate(contrast = "aPS_vs_PD"))
  topn_terms = topn_terms %>% rownames_to_column("terms")
  write.csv(topn_terms,quote = F,row.names = F, file = file.path(output_dir,paste0(sub_name, "_ssGSEA_topn_terms_for_heatmap.csv")))
  
  
  # Order the columns of 'plot_df' according to 'annolabel'
  annolabel <- data.frame(Group = group_info[,"Condition"])
  rownames(annolabel) <- group_info$ID
  # Order the columns of 'plot_df' according to 'annolabel'
  ordered_columns <- rownames(annolabel)[order(annolabel$Condition)]
  group_colors <- list(Condition = c("aPS" = "#8DA0CB","Ctrl" = "#FC8D62", "PD" = "#A6D854"))
  
  plot_df <- gsva_results[unique(c(rownames(aPS_vs_Ctrl),rownames(PD_vs_Ctrl),rownames(aPS_vs_PD))), ordered_columns]
  
  #legend_breaks <- seq(from = round(min(plot_df), 1), to = round(max(plot_df), 1), by = 0.2)
  
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
  dev.off()
  
  ### Means heatmap
  Condition_means <- gsva_results %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    mutate(Condition = annolabel$Condition) %>%
    group_by(Condition) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    column_to_rownames(var = "Condition") %>%
    t()
  
  write_Condition_means = Condition_means %>% as.data.frame() %>%rownames_to_column("terms")
  write.csv(write_Condition_means,quote = F,row.names = F, file = file.path(output_dir, paste0(sub_name,"_ssGSEA_Condition_means.csv")))
  
  Condition_means = Condition_means[unique(c(rownames(aPS_vs_Ctrl),rownames(PD_vs_Ctrl),rownames(aPS_vs_PD))),]
  
  # Custom function to abbreviate long terms by adding "..." at the end
  rownames(Condition_means) <- sapply(rownames(Condition_means), abbreviate_terms_end)
  pathway_heatmap_mean <- pheatmap::pheatmap(Condition_means,
                                             fontsize_row = 6, cluster_cols = FALSE, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), border_color = NA)
  
  
  pdf(file.path(output_dir,paste0(sub_name,"_ssGSEA_Condition_means.pdf")),width = 7)
  print(pathway_heatmap_mean)
  dev.off()
  
}
