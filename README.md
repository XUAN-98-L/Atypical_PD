Bioinformatics analysis of PD CSF proteomics

Based on the protein expression matrix and ontology gene sets from the Molecular Signature Database (MSigDB v7.5.1) with at least 10 overlapping genes, the single-sample gene set enrichment analysis (ssGSEA) algorithm was used to perform gene ontology analysis by using R package “gsva”.

The ssGSEA scores were normalized using Z-scores and subsequently used as input for the R package limma to estimate the significance of pathways among subgroups through moderated t-statistics. Pathways with Benjamini-Hochberg corrected p-values below 0.1 were considered significantly enriched. The top 20 pathways, ranked by t-statistic for each comparison group, were visualized using a heatmap.

The terms in ssGSEA_topn_terms_for_heatmap.csv are visualized in heatmaps.
The file ssGSEA_topn_terms_*_with_matched_proteins.xlsx provides the proteins matched to each geneset.

Comparison group:
Results_NPS_Ctr_PD:
NPS_vs_PD
NPS_vs_Ctr

Results_aPS_Ctrl_PD:
aPS_vs_Ctrl
PD_vs_Ctrl

