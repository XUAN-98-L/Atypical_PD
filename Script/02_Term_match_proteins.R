
library(dplyr)
library(tidyverse)
library(biomaRt)
# Loading data
expression_data = readxl::read_excel("Data/CSF_APS_individual_patient_results.xlsx",sheet = 1)
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

metadata = readxl::read_excel("Data/CSF_APS_individual_patient_results.xlsx",sheet = 2)
group_info = metadata %>% dplyr::select("ID","Subtype")
# Focus on subtype excluded Tauopathy
group_info = group_info %>% dplyr::filter(Subtype != "Tauopathy")
# replace "Non-PD-Synucleinopathy" into "NPS"
group_info$Subtype[group_info$Subtype == "Non-PD-Synucleinopathy"] = "NPS"

protein_expression = protein_expression[,group_info$ID]


GO_gene_sets <- msigdbr::msigdbr(
    species = "Homo sapiens", # Can change this to what species you need
    category = "C5" # Only GO gene sets
  ) 
GO_gene_sets = GO_gene_sets[which(GO_gene_sets$gs_subcat =="GO:BP"),]
# Background removal
GO_gene_sets <- GO_gene_sets[GO_gene_sets$gene_symbol %in% rownames(protein_expression), ]

# Split geneset
GO_list <- split(
  GO_gene_sets$gene_symbol, # The genes we want split into pathways
  GO_gene_sets$gs_name # The pathways made as the higher levels of the list
)

# Extract gene names from expression matrix
genes_of_interest <- rownames(protein_expression)

# Function to calculate overlap with genes of interest
overlap_with_genes <- function(geneset) {
  length(intersect(geneset, genes_of_interest))
}
# Filter GO_list, keeping only sets with at least 10 overlapping genes
filtered_GO_list <- GO_list[sapply(GO_list, overlap_with_genes) >= 10]

# Convert filtered_terms to a data frame
filtered_df <- data.frame(
  Term = names(filtered_GO_list),
  Genes = sapply(filtered_GO_list, function(genes) paste(unique(genes), collapse = ","))
)

output_dir = "Results/02_Term_match_proteins"
if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}
write.csv(filtered_df, file = file.path(output_dir, "Term_match_proteins.csv"), row.names = FALSE)