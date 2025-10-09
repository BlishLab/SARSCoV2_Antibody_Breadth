

#Load libraries
library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(MoMAColors)
message("finished libraries")

#Clear current environment (to ensure code runs only with what is given here)
rm(list = ls())

#Set data path and figure path
data_path <- ""

fig_path <- ""
message("finished data/fig paths")

# Convenience functions
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}
message("finished convenience functions")

#Load ligand-receptor network and ligand-target matrix
organism = "human"
lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
ligand_target_matrix = readRDS("/oak/stanford/groups/cblish/Rebecca/Parse/analysis/IRIS_megakit_analysis_code/ligand_target_matrix_nsga2r_final.rds")
colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
message("finished line85")

#Load seurat object, set idents for research question
seu <- readRDS("scPBMCcovid.rds")
Idents(seu) <- "Breadth"

#convert to SCE, update gene symbols
sce <- Seurat::as.SingleCellExperiment(seu, assay = "RNA")
sce <- alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()
message("finished line94")

#Define metadata columns
SummarizedExperiment::colData(sce)$Donor = SummarizedExperiment::colData(sce)$Donor %>% make.names()
SummarizedExperiment::colData(sce)$manual.coarse.idents = SummarizedExperiment::colData(sce)$manual.coarse.idents %>% make.names()
SummarizedExperiment::colData(sce)$Breadth = SummarizedExperiment::colData(sce)$Breadth %>% make.names()
message("madenames")
sample_id = "Donor"
group_id = "Breadth"
celltype_id = "manual.coarse.idents"
covariates = NA
batches = NA

#Define sender and receiver cell types
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = c("NK")

message("finished line109")

#Extract cell type abundance, expression data
min_cells = 10
abundance_expression_info = get_abundance_expression_info(sce = sce, sample_id = sample_id, group_id = group_id,
                                                          celltype_id = celltype_id, min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi,
                                                          lr_network = lr_network)
abundance_plot <- abundance_expression_info$abund_plot_sample
abundance_plot
SaveFigure(abundance_plot, "Breadth_abundance_plot", width = 9, height = 18)
message("finished saving abundance plot")


#Define contrasts, covariates, perform DE analysis
contrasts_oi = c("'Narrow-Broad','Broad-Narrow'")
contrast_tbl = tibble(
  contrast = c("Narrow-Broad","Broad-Narrow"),
  group = c("Narrow","Broad"))
DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches,
                      covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells)
DE_info$celltype_de$de_output_tidy %>% arrange(p_adj) %>% head()
celltype_de = DE_info$celltype_de$de_output_tidy
message("finished DE")

#Combine DE for ligand-senders and receptors-receivers
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)
sender_receiver_de %>% head(20)
message("finished senderrecieverDE")

#Predict NicheNet ligand activities and ligand-target links based on DE results
logFC_threshold = 0.25
p_val_threshold = 0.05
fraction_cutoff = 0.05
p_val_adj = F #In case of more samples per group, and high number of DE genes per group (>50), recommend using adjusted p values
top_n_target = 250 #select top n of predicted target genes to be considered
verbose = TRUE
cores_system = 8
n.cores = min(cores_system, sender_receiver_de$receiver %>% unique() %>% length()) # use one core per receiver cell type
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target,
  verbose = verbose,
  n.cores = n.cores
)))
message("finished activities")

#Check DE genes used for activity analysis and the output of the analysis
ligand_activities_targets_DEgenes$de_genes_df %>% head(20)
ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

#Define prioritization weights, prepare grouping objects
prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)
prioritizing_weights = c(prioritizing_weights_DE,
                         prioritizing_weights_activity,
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency,
                         prioritizing_weights_relative_abundance)
sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()
message("finished prioritization")

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}
message("finished line199")

#Run prioritization
prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  prioritizing_weights = prioritizing_weights,
  fraction_cutoff = fraction_cutoff,
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender
))
message("finished line214")

prioritization_tables$group_prioritization_tbl %>% head(20)

lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(),
                                                    abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix,
                                                    logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj)
message("finished line221")

#Save output of MultiNicheNet
path = ""

multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
)
multinichenet_output = make_lite_output(multinichenet_output)
message("finished line237")

save = TRUE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output_NKAll.rds"))
  
}
message("finished line244 - saved")

