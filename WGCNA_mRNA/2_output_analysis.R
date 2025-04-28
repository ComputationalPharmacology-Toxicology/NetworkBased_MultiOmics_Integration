## SECOND PART OF WGCNA: analysis of the results ##

# Load libraries
library(devtools)
library(WGCNA)
library(flashClust)
library(curl)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(CorLevelPlot)
library(pheatmap)
library(data.table)


## 0. Definition of input/output paths ----------------------------------------------------------------------
initial_path <- file.path("WGCNA_mRNA/data")
out_path <- file.path("WGCNA_mRNA/output")
plot_path <- file.path("WGCNA_mRNA/plot")



## 1. Load and explor output ----------------------------------------------------------------------
bwnet <- readRDS(file.path(out_path, "bwnet.rds"))
str(bwnet)
# Module Eigengenes 
module_eigengenes <- bwnet$MEs
# Print out a preview
head(module_eigengenes)

# Get number of genes for each module
num.per.module <- table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
pdf(file.path(plot_path, '2_dendro_colors_plot.pdf'))
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)
dev.off()


## 2A. Relate modules to traits (=treatments) --------------------------------------------------

# Create colData with treatments
sampletype <- as.factor(c(rep('CTR',6), rep('GEMTAX', 5), rep('SIM',6), rep('VPA', 6), 
rep('VPA_SIM', 5), rep('VS_GEMTAX', 6)))
colData <- data.frame(sampletype, row.names = colnames(data))

# Create traits file - binarize categorical variables: ALL PAIRWISE COMPARISONS
traits_tot <- binarizeCategoricalColumns(colData$sampletype,
                           includePairwise = TRUE,
                           includeLevelVsAll = TRUE,
                           dropFirstLevelVsAll = FALSE,
                           minCount = 1)
rownames(traits_tot) <- rownames(module_eigengenes)

# Define numbers of genes and samples
norm.counts <- readRDS(file.path(out_path, "norm_counts.rds"))
nSample_tot <- nrow(norm.counts)
nGenes_tot <- ncol(norm.counts)

# Compute correlation between eigengenes of modules and traits (binarized treatments)
module.trait_tot.corr <- cor(module_eigengenes, traits_tot, use = 'p')
module.trait_tot.corr.pvals <- corPvalueStudent(module.trait_tot.corr, nSample_tot)

# Visualize module-trait association as a heatmap
heatmap.data_tot <- merge(module_eigengenes, traits_tot, by = 'row.names')
heatmap.data_tot <- heatmap.data_tot %>% 
  column_to_rownames(var = 'Row.names')
head(heatmap.data_tot)

# Change column names by adding the number of genes per module
colnames(heatmap.data_tot) <- sapply(colnames(heatmap.data_tot), function(col) {
  module_name_tot <- gsub("^ME", "", col)
  if (module_name_tot %in% names(num.per.module)) {
    return(paste0(col, " (", num.per.module[module_name_tot], ")"))
  } else {
    return(col) 
  }
})
 
colnames(heatmap.data_tot)

# Create the heatmap with all the pairwise combinations
pdf(file.path(plot_path, '2_heatmap_modules_vs_traits_tot.pdf'), width = 16, height = 12)
CorLevelPlot(heatmap.data_tot,
             x = names(heatmap.data_tot)[grep("^data", names(heatmap.data_tot))],
             y = names(heatmap.data_tot)[grep("^ME", names(heatmap.data_tot))],
             col = c("blue1", "skyblue", "white", "pink", "red"),
             signifCutpoints = c(0, 0.001, 0.01, 0.03, 1), 
             rotLabX = 90)
dev.off()

### From now on: KEEP ONLY TRAIT VS ALL ###
# Create traits file - binarize categorical variables
traits <- binarizeCategoricalColumns(colData$sampletype,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE,
                           dropFirstLevelVsAll = FALSE,
                           minCount = 1)
rownames(traits) <- rownames(module_eigengenes)

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

# Compute correlation between eigengenes of modules and traits (binarized treatments)
module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# Visualize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')
head(heatmap.data)

# Change column names by adding the number of genes per module
colnames(heatmap.data) <- sapply(colnames(heatmap.data), function(col) {
  module_name <- gsub("^ME", "", col)
  if (module_name %in% names(num.per.module)) {
    return(paste0(col, " (", num.per.module[module_name], ")"))
  } else {
    return(col) 
  }
})
 
# Create the heatmap with only "vs all"
colnames(heatmap.data)

pdf(file.path(plot_path, '2_heatmap_modules_vs_traits.pdf'))
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[grep("^data", names(heatmap.data))],
             y = names(heatmap.data)[grep("^ME", names(heatmap.data))],
             col = c("blue1", "skyblue", "white", "pink", "red"),
             signifCutpoints = c(0, 0.001, 0.01, 0.03, 1), 
             rotLabX = 90)
dev.off()

# Extract genes for modules
module.gene.mapping <- as.data.frame(bwnet$colors)

# Modules selection: keep only row with p-value < pval_star_threshold
pval_star_threshold <- 0.01
significant_rows <- rownames(module.trait.corr.pvals)[apply(module.trait.corr.pvals < pval_star_threshold, 1, any)]
colors <- sub("^ME", "", significant_rows)
colors

saveRDS(colors, file = file.path(out_path, "colors.rds"))

# For each significant row, keep only columns with p-value < pval_star_threshold
significant_cols <- lapply(significant_rows, function(row) {
  which(module.trait.corr.pvals[row, ] < pval_star_threshold)
})
names(significant_cols) <- colors
significant_cols

# Load txt file with the number of reads
file = file.path(initial_path, 'salmon.merged.gene_counts.tsv')
# Read the file and skip the first line and read only the columns "gene_id" e "gene_name"
counts_data <- fread(file, sep = "\t", select = c("gene_id", "gene_name"))
head(counts_data)

# Merge of module.gene.mapping to obtain the gene symbols
head(module.gene.mapping)
module.gene.mapping$gene_id <- rownames(module.gene.mapping)
merged_data <- merge(module.gene.mapping, counts_data, by = "gene_id", all.x = TRUE)
dim(module.gene.mapping)[1] == dim(merged_data)[1]
# Check NAs in gene_name
print(sum(is.na(merged_data$gene_name)))
# Create combined name
merged_data$combined_name <- paste0(merged_data$gene_name, "_", merged_data$gene_id)


# Create and save gene lists
gene_lists <- list()
gene_lists_comb_names <- list()
for (color in colors) {
  gene_lists[[color]] <- merged_data %>% 
    filter(`bwnet$colors` == color) %>%
    pull(gene_id)
  gene_lists_comb_names[[color]] <- merged_data %>% 
    filter(`bwnet$colors` == color) %>%
    pull(combined_name)
  print(dim(unique(gene_lists[[color]]))==dim(gene_lists[[color]]))
  writeLines(unique(gene_lists[[color]]), file.path(out_path, paste0(color, "_genes.txt")))
}

# Print dimentions of modules
lapply(gene_lists, length)


## 2B. Intramodular analysis: Identifying driver genes --------------------------------------------------

# Create cobined name to have gene symbols in the plots later
counts_data$combined_name <- paste0(counts_data$gene_name, "_", counts_data$gene_id)
t_norm.counts <- t(norm.counts)
t_norm_df <- as.data.frame(t_norm.counts) #from matrix to dataframe
t_norm_df$gene_id <- rownames(t_norm.counts)
# Merge to get combined_name
merged_df <- merge(t_norm_df, counts_data[, c("gene_id", "combined_name")], 
                   by = "gene_id", all.x = TRUE)
# Use combined_name as new rownames
rownames(merged_df) <- merged_df$combined_name
# Remove extra columns
merged_df$gene_id <- NULL
merged_df$combined_name <- NULL
norm.counts <- t(as.matrix(merged_df)) 

# Compute correlation between normalized counts and traits (and compute adjusted pvalues)
# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.
module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)
# See results
module.membership.measure[,1:10]
module.membership.measure.pvals[,1:10] 


# Compute correlation between normalized counts and traits (and compute adjusted pvalues)
# Define a list with gene set names and conditions to be used for filtering
gene_groups <- lapply(names(significant_cols), function(color) {
  list(genes = gene_lists[[color]], conditions = significant_cols[[color]])
})
names(gene_groups) <- names(significant_cols)

# Repeat the same but with combined_names
gene_groups_comb_names <- lapply(names(significant_cols), function(color) {
  list(genes = gene_lists_comb_names[[color]], conditions = significant_cols[[color]])
})
names(gene_groups_comb_names) <- names(significant_cols)


# i) Create heatmap for each ME module 
# Here we consider genes significant for at least one of the selected traits
for (group_name in names(gene_groups_comb_names)) {
  
  # Get genes and conditions
  genes <- gene_groups_comb_names[[group_name]]$genes
  conditions <- gene_groups_comb_names[[group_name]]$conditions
  
  # Calculate the correlation matrix
  gene.signf.corr <- cor(norm.counts[, genes], traits, use = 'p')
  
  # Calcultate p-value matrix with corPvalueStudent function 
  gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
  # column-wise p-value correction 
  gene.signf.corr.padj <- apply(gene.signf.corr.pvals, 2, function(x) p.adjust(x, method = "BH"))

  # Selection of significant genes
  p_threshold <- 0.01
  
  if (length(conditions) == 1) {
    selected_genes <- which(gene.signf.corr.padj[, conditions] < p_threshold)
  } else {
    selected_genes <- which(rowSums(gene.signf.corr.padj[, conditions] < p_threshold) > 0)
  }

  cat(paste0('\n', group_name, ': ', length(selected_genes)))

  # Check if there are selected genes
  if (length(selected_genes) > 0) {
    # Filter correlation and annotation matrix
    filtered_corr <- gene.signf.corr[selected_genes, , drop = FALSE]
    filtered_padj <- gene.signf.corr.padj[selected_genes, , drop = FALSE]
    
    # Create annotations for the heatmap
    annot_matrix <- matrix(paste0(round(filtered_corr, 2), " (", signif(filtered_padj, 2), ")"), 
                           nrow = nrow(filtered_corr), ncol = ncol(filtered_corr))
    rownames(annot_matrix) <- rownames(filtered_corr)
    colnames(annot_matrix) <- colnames(filtered_corr)
    
    # Symmetrical color scale between -max_corr and max_corr
    color_scale <- colorRampPalette(c("blue", "white", "red"))(50)
    
    # Save heatmap with only filtered genes
    pdf(file.path(plot_path, paste0('2_heatmap_', group_name, '.pdf')), 
    height = max(5, nrow(filtered_corr) * 0.2), width = 10)  
    pheatmap(filtered_corr, 
             display_numbers = annot_matrix,  
             color = color_scale,  
             cluster_rows = FALSE, cluster_cols = FALSE, 
             fontsize_number = 10,  
             cellwidth = 80,  
             fontface = 'bold',
             number_color = 'black',
             main = paste0("Heatmap Correlazione (p.adj < ", p_threshold, " in ", group_name, ")"),
             breaks = seq(-1, 1, length.out = 51))  
    dev.off()
  }
}


# ii) Create txt files with gene names for each ME module and each significant condition  
for (group_name in names(gene_groups_comb_names)) {
  
  # Take genes and conditions
  genes <- gene_groups_comb_names[[group_name]]$genes
  conditions <- gene_groups_comb_names[[group_name]]$conditions
  
  # Calculate the correlation matrix
  gene.signf.corr <- cor(norm.counts[, genes], traits, use = 'p')
  # Calculate the p-value matrix with corPvalueStudent function
  gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
  # column-wise p-value correction 
  gene.signf.corr.padj <- apply(gene.signf.corr.pvals, 2, function(x) p.adjust(x, method = "BH"))
  p_threshold <- 0.01
  
  # Selection of significant genes
  for (i in seq_along(conditions)) {
    condition <- conditions[[i]] 
    name <- names(conditions)[i]  
    
    selected_genes <- rownames(gene.signf.corr.padj)[which(gene.signf.corr.padj[, condition] < p_threshold)]
    selected_genes_clean <- sub(".*_", "", selected_genes)
    
    # Save the list in a txt file 
    if (length(selected_genes_clean) > 0) {
      cat(paste0('\n', group_name, ' - ', name, ': ', length(selected_genes_clean)))
      file_name <- file.path(out_path, paste0("significant_genes_", group_name, "_", name, ".txt"))
      writeLines(selected_genes_clean, file_name)
    }
  }
}
