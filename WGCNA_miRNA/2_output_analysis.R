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
library(RColorBrewer)
library(readxl)
library(dplyr)
library(purrr)


## 0. Definition of input/output paths ----------------------------------------------------------------------
initial_path <- file.path("WGCNA_miRNA/data")
out_path <- file.path("WGCNA_miRNA/output")
plot_path <- file.path("WGCNA_miRNA/plot")



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
rep('VPASIM', 5), rep('VSGT', 6)))
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
nSamples_tot <- nrow(norm.counts)
nGenes_tot <- ncol(norm.counts)

# Compute correlation between eigengenes of modules and traits (binarized treatments)
module.trait_tot.corr <- cor(module_eigengenes, traits_tot, use = 'p')
module.trait_tot.corr.pvals <- corPvalueStudent(module.trait_tot.corr, nSamples_tot)

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
                           includePairwise = FALSE, # not all pairwise comparisons
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

colnames(heatmap.data)

# Create the heatmap with only "vs all"
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
pval_star_threshold <- 0.03
significant_rows <- rownames(module.trait.corr.pvals)[apply(module.trait.corr.pvals < pval_star_threshold, 1, any)]
colors <- sub("^ME", "", significant_rows)# remove "ME" from names
colors

saveRDS(colors, file = file.path(out_path, "colors.rds"))

# For each significant row, keep only columns with p-value < pval_star_threshold
significant_cols <- lapply(significant_rows, function(row) {
  which(module.trait.corr.pvals[row, ] < pval_star_threshold)
})
names(significant_cols) <- colors
significant_cols <- significant_cols[names(significant_cols) != "grey"]
significant_cols

# Note: we decided to exclude the gray module since having other modules more relevant to that trait
colors <- colors[colors != "grey"]
saveRDS(colors, file = file.path(out_path, "colors.rds"))

# Create and save gene lists
gene_lists <- list()
for (color in colors) {
  gene_lists[[color]] <- module.gene.mapping %>% 
    filter(`bwnet$colors` == color) %>%
    rownames()
  writeLines(gene_lists[[color]], file.path(out_path, paste0(color, "_genes.txt")))
}

# Print dimentions of modules
lapply(gene_lists, length)



## 2B. Intramodular analysis: Identifying driver genes --------------------------------------------------

# Compute correlation between eigengenes of modules and normalized counts
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


# i) Create heatmap for each ME module 
# Here we consider genes significant for at least one of the selected traits
for (group_name in names(gene_groups)) {

  # Get genes and conditions
  genes <- gene_groups[[group_name]]$genes
  conditions <- gene_groups[[group_name]]$conditions
  
  # Calculate the correlation matrix
  gene.signf.corr <- cor(norm.counts[, genes], traits, use = 'p')
  
  # Calcultate p-value matrix with corPvalueStudent function 
  gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
  # column-wise p-value correction 
  gene.signf.corr.padj <- apply(gene.signf.corr.pvals, 2, function(x) p.adjust(x, method = "BH"))
  
  # Selection of significant genes
  p_threshold <- 0.05
  
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

# ii) Heatmaps of norm counts for driver genes in each MEs 
# Load txt file with the number of reads
file = file.path(initial_path, 'mature_counts.csv')

# Read the file and skip the first line
counts_data <- read.csv(file, sep = ",", header = TRUE, row.names=1)
counts_data <- counts_data[order(rownames(counts_data)), ] # order samples by name
rownames(counts_data)
counts_data <- t(counts_data) # get the transposed matrix of counts_data
head(counts_data)

norm_counts <- log10(counts_data)
norm_counts[norm_counts == -Inf] <- -1000
head(norm_counts)

rownames(colData) <- rownames(module_eigengenes)
heat_colors <- brewer.pal(6, "YlOrRd")


# Create txt files and a heatmap plot with norm counts (norm: log10) with gene names for each ME module and each significant condition 
for (group_name in names(gene_groups)) {
  
  # Take genes and conditions
  genes <- gene_groups[[group_name]]$genes
  conditions <- gene_groups[[group_name]]$conditions
  
  # Calculate the correlation matrix
  gene.signf.corr <- cor(norm.counts[, genes], traits, use = 'p')
  # Calculate the p-value matrix with corPvalueStudent function
  gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
  # column-wise p-value correction 
  gene.signf.corr.padj <- apply(gene.signf.corr.pvals, 2, function(x) p.adjust(x, method = "BH"))
  p_threshold <- 0.05
  
  # Selection of significant genes
  for (i in seq_along(conditions)) {
    condition <- conditions[[i]]  
    name <- names(conditions)[i]  
    
    selected_genes <- rownames(gene.signf.corr.padj)[which(gene.signf.corr.padj[, condition] < p_threshold)]
    
    # Save the list in a txt file and create the heatmap
    if (length(selected_genes) > 1) { 
      cat(paste0('\n', group_name, ' - ', name, ': ', length(selected_genes)))
      file_name <- file.path(out_path, paste0("significant_genes_", group_name, "_", name, ".txt"))
      writeLines(selected_genes, file_name)

      # Save the norm counts heatmaps only for significant genes
      norm_counts_color <- norm_counts[selected_genes, ]
      pdf(file.path(plot_path, paste0('2_heatmap_norm_counts_', group_name, '_', name, '.pdf')), height = max(5, nrow(norm_counts_color) * 0.2), width = 15)
      pheatmap(norm_counts_color, 
            color = heat_colors, 
            cluster_rows = T, 
            cluster_cols = F, 
            show_rownames = T,
            annotation = colData, 
            border_color = NA, 
            fontsize = 10, 
            scale = "row", 
            fontsize_row = 8, 
            cellwidth = 20,
            main = paste0("Heatmap Norm Counts - ", group_name))
      dev.off()
    }
  }
}
