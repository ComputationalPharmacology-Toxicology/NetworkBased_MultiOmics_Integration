## SECOND PART OF WGCNA: analysi of the results

# create plots and graph

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

# # load output
bwnet <- readRDS(file.path(out_path, "bwnet.rds"))
str(bwnet)

# 2. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs
# Print out a preview
head(module_eigengenes)

# get number of genes for each module
num.per.module <- table(bwnet$colors)
    # black      blue     brown     green      grey   magenta      pink    purple 
    #    61       131       124        83       593        29        41        23 
    #   red turquoise    yellow 
    #    77       158        99 

# Plot the dendrogram and the module colors before and after merging underneath
pdf(file.path(plot_path, '2_dendro_colors_plot.pdf'))
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)
dev.off()


# 3A. Relate modules to traits --------------------------------------------------

# Create colData
sampletype <- as.factor(c(rep('CTR',6), rep('GEMTAX', 5), rep('SIM',6), rep('VPA', 6), 
rep('VPASIM', 5), rep('VSGT', 6)))
colData <- data.frame(sampletype, row.names = colnames(data))
# create traits file - binarize categorical variables: ALL PAIRWISE COMPARISONS
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

# qui lo calcoliamo ma in realtà la funzione della heatmap calcola la correlazione e i pvalue direttamente
# compute correlation between eigengenes of modules and traits (binarized treatments)
module.trait_tot.corr <- cor(module_eigengenes, traits_tot, use = 'p')
module.trait_tot.corr.pvals <- corPvalueStudent(module.trait_tot.corr, nSamples_tot)

# visualize module-trait association as a heatmap
heatmap.data_tot <- merge(module_eigengenes, traits_tot, by = 'row.names')
heatmap.data_tot <- heatmap.data_tot %>% 
  column_to_rownames(var = 'Row.names')
head(heatmap.data_tot)

# Modifica i nomi delle colonne
colnames(heatmap.data_tot) <- sapply(colnames(heatmap.data_tot), function(col) {
  # Cerca il nome del modulo rimuovendo "ME"
  module_name_tot <- gsub("^ME", "", col)
  # Se il modulo è nella lista, aggiungi la numerosità
  if (module_name_tot %in% names(num.per.module)) {
    return(paste0(col, " (", num.per.module[module_name_tot], ")"))
  } else {
    return(col)  # Se non è un modulo, lascia il nome invariato
  }
})
 
# Controlla i nuovi nomi delle colonne
colnames(heatmap.data_tot)

pdf(file.path(plot_path, '2_heatmap_modules_vs_traits_tot.pdf'), width = 16, height = 12)
CorLevelPlot(heatmap.data_tot,
             x = names(heatmap.data_tot)[grep("^data", names(heatmap.data_tot))],
             y = names(heatmap.data_tot)[grep("^ME", names(heatmap.data_tot))],
             col = c("blue1", "skyblue", "white", "pink", "red"),
             signifCutpoints = c(0, 0.001, 0.01, 0.03, 1), 
             rotLabX = 90)
dev.off()

### KEEP ONLY TRAIT VS ALL ###
# create traits file - binarize categorical variables
traits <- binarizeCategoricalColumns(colData$sampletype,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE,
                           dropFirstLevelVsAll = FALSE,
                           minCount = 1)
rownames(traits) <- rownames(module_eigengenes)

# Define numbers of genes and samples
# norm.counts <- readRDS(file.path(out_path, "norm_counts.rds"))
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

# qui lo calcoliamo ma in realtà la funzione della heatmap calcola la correlazione e i pvalue direttamente
# compute correlation between eigengenes of modules and traits (binarized treatments)
module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')
head(heatmap.data)

# Modifica i nomi delle colonne
colnames(heatmap.data) <- sapply(colnames(heatmap.data), function(col) {
  # Cerca il nome del modulo rimuovendo "ME"
  module_name <- gsub("^ME", "", col)
  # Se il modulo è nella lista, aggiungi la numerosità
  if (module_name %in% names(num.per.module)) {
    return(paste0(col, " (", num.per.module[module_name], ")"))
  } else {
    return(col)  # Se non è un modulo, lascia il nome invariato
  }
})
 
# Controlla i nuovi nomi delle colonne
colnames(heatmap.data)

pdf(file.path(plot_path, '2_heatmap_modules_vs_traits.pdf'))
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[grep("^data", names(heatmap.data))],
             y = names(heatmap.data)[grep("^ME", names(heatmap.data))],
             col = c("blue1", "skyblue", "white", "pink", "red"),
             signifCutpoints = c(0, 0.001, 0.01, 0.03, 1), 
             rotLabX = 90)
dev.off()

# extract genes for modules
module.gene.mapping <- as.data.frame(bwnet$colors)

# selezione moduli
# Seleziona le righe in cui almeno un p-value è < pval_star_threshold 
# NB qui il modulo grey resta nella lista ma non lo guardiamo nelle analisi successive perchè
# per quel confronto abbiamo già 4 moduli molto più correlati
# se lo escludessimo qui con la soglia dovremmo escludere anche il modulo black
pval_star_threshold <- 0.03
significant_rows <- rownames(module.trait.corr.pvals)[apply(module.trait.corr.pvals < pval_star_threshold, 1, any)]
# Rimuove il prefisso "ME" dai nomi
colors <- sub("^ME", "", significant_rows)
colors


# Per ogni riga significativa, trova le colonne con p-value < pval_star_threshold
significant_cols <- lapply(significant_rows, function(row) {
  which(module.trait.corr.pvals[row, ] < pval_star_threshold)
})
names(significant_cols) <- colors
significant_cols <- significant_cols[names(significant_cols) != "grey"]
significant_cols

# escludiamo gray perchè non siamo interessati essendo poco significativo e avendo altri moduli più rilevanti per quel tratto
colors <- colors[colors != "grey"]
saveRDS(colors, file = file.path(out_path, "colors.rds"))

# Creazione e salvataggio delle liste di geni in un unico ciclo
gene_lists <- list()
for (color in colors) {
  gene_lists[[color]] <- module.gene.mapping %>% 
    filter(`bwnet$colors` == color) %>%
    rownames()
  # Salvataggio immediato
  writeLines(gene_lists[[color]], file.path(out_path, paste0(color, "_genes.txt")))
}

# print dimentions of modules
lapply(gene_lists, length)

# $black
# [1] 61

# $yellow
# [1] 99

# $red
# [1] 77

# $blue
# [1] 131

# $brown
# [1] 124


# 3B. Intramodular analysis: Identifying driver genes ---------------

# compute correlation between eigengenes of modules and normalized counts

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.
module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure[,1:10]
module.membership.measure.pvals[,1:10] 


# compute correlation between normalized counts and traits (and compute adjusted pvalues)
#( nel nostro caso penso che possiamo farlo per ognuna delle condizioni, quindi vedere quali sono i 
# geni più associati al trattamento con sim con valp o entrambi)

# Definisci una lista con i nomi dei gruppi di geni e le condizioni da utilizzare per il filtro
gene_groups <- lapply(names(significant_cols), function(color) {
  list(genes = gene_lists[[color]], conditions = significant_cols[[color]])
})
# Assegna i nomi ai gruppi di geni
names(gene_groups) <- names(significant_cols)

### Create heatmap for each ME module ###
#(qui consideriamo i geni significativi per almeno uno dei tratti selezionati)
for (group_name in names(gene_groups)) {
  
  # Prendi i geni e le condizioni
  genes <- gene_groups[[group_name]]$genes
  conditions <- gene_groups[[group_name]]$conditions
  
  # Calcola la matrice di correlazione
  gene.signf.corr <- cor(norm.counts[, genes], traits, use = 'p')
  
  # Calcola la matrice di p-value con corPvalueStudent
  gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
  
  # Correggi i p-value per colonna mantenendo la struttura originale
  gene.signf.corr.padj <- apply(gene.signf.corr.pvals, 2, function(x) p.adjust(x, method = "BH"))

  p_threshold <- 0.05
  # Seleziona i geni con padj < p_threshold nelle condizioni specificate
  if (length(conditions) == 1) {
    # Se c'è una sola condizione, trasformiamo in una matrice
    selected_genes <- which(gene.signf.corr.padj[, conditions] < p_threshold)
  } else {
    # Se ci sono più condizioni, usa rowSums
    selected_genes <- which(rowSums(gene.signf.corr.padj[, conditions] < p_threshold) > 0)
  }

  cat(paste0('\n', group_name, ': ', length(selected_genes)))

  # Verifica se ci sono geni selezionati
  if (length(selected_genes) > 0) {
    # Filtra la matrice di correlazione e annotazioni
    filtered_corr <- gene.signf.corr[selected_genes, , drop = FALSE]
    filtered_padj <- gene.signf.corr.padj[selected_genes, , drop = FALSE]
    
    # Crea annotazioni per la heatmap
    annot_matrix <- matrix(paste0(round(filtered_corr, 2), " (", signif(filtered_padj, 2), ")"), 
                           nrow = nrow(filtered_corr), ncol = ncol(filtered_corr))
    rownames(annot_matrix) <- rownames(filtered_corr)
    colnames(annot_matrix) <- colnames(filtered_corr)
    
    # Scala di colori simmetrica tra -max_corr e max_corr
    color_scale <- colorRampPalette(c("blue", "white", "red"))(50)
    
    # Salva la heatmap solo con i geni filtrati
    pdf(file.path(plot_path, paste0('2_heatmap_', group_name, '.pdf')), 
    height = max(5, nrow(filtered_corr) * 0.2), width = 10)  # Altezza dinamica
    pheatmap(filtered_corr, 
             display_numbers = annot_matrix,  # Mostra solo i geni filtrati
             color = color_scale,  # Scala di colori simmetrica
             cluster_rows = FALSE, cluster_cols = FALSE, # Clustering automatico
             fontsize_number = 10,  # Dimensione del testo nelle celle
             cellwidth = 80,  # Larghezza delle celle
             fontface = 'bold',
             number_color = 'black',
             main = paste0("Heatmap Correlazione (p.adj < ", p_threshold, " in ", group_name, ")"),
             breaks = seq(-1, 1, length.out = 51))  # Imposta i limiti della scala di colori
    dev.off()
  }
}

########## Heatmaps of norm counts for driver genes in each MEs ##########
# Load txt file created by nf-core pipeline, with the number of reads for each mirna 
file = file.path(initial_path, 'mature_counts.csv')
# Read the file and skip the first line
counts_data <- read.csv(file, sep = ",", header = TRUE, row.names=1)#, colClasses = c(rep("character",2), rep("integer", 12)))
counts_data <- counts_data[order(rownames(counts_data)), ] # order samples by name
rownames(counts_data)
counts_data <- t(counts_data) # Get the transposed matrix of counts_data
head(counts_data)

norm_counts <- log10(counts_data)
norm_counts[norm_counts == -Inf] <- -1000
head(norm_counts)

## Options for pheatmap
# Annotate heatmap (optional)
rownames(colData) <- rownames(module_eigengenes)
# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")


### Create txt files with gene names for each ME module and each significant condition ###
#(if a module is sign for two different condition, we will have 2 txt) 
### Create heatmap plots with norm counts (log10) with gene names for each ME module and each significant condition ###
for (group_name in names(gene_groups)) {
  
  # Prendi i geni e le condizioni
  genes <- gene_groups[[group_name]]$genes
  conditions <- gene_groups[[group_name]]$conditions
  
  # Calcola la matrice di correlazione
  gene.signf.corr <- cor(norm.counts[, genes], traits, use = 'p')
  # Calcola la matrice di p-value con corPvalueStudent
  gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
  # Correggi i p-value per colonna mantenendo la struttura originale
  gene.signf.corr.padj <- apply(gene.signf.corr.pvals, 2, function(x) p.adjust(x, method = "BH"))
  p_threshold <- 0.05
  
  # Seleziona i geni con padj < p_threshold nelle condizioni specificate
  for (i in seq_along(conditions)) {
    condition <- conditions[[i]]  # Ottieni la condizione
    name <- names(conditions)[i]  # Ottieni il nome della condizione
    
    selected_genes <- rownames(gene.signf.corr.padj)[which(gene.signf.corr.padj[, condition] < p_threshold)]
    
    # Salva la lista in un file esterno e fai heatmap solo se la lista contiene geni
    if (length(selected_genes) > 1) { # maggiore stretto di 1 perchè altrimenti pheatmap crea problemi essendo solo un vettore (numeric) e non una matrice
      # salva lista geni
      cat(paste0('\n', group_name, ' - ', name, ': ', length(selected_genes)))
      file_name <- file.path(out_path, paste0("significant_genes_", group_name, "_", name, ".txt"))
      writeLines(selected_genes, file_name)

      # Salva la heatmap per norm counts solo con i geni significativi
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
