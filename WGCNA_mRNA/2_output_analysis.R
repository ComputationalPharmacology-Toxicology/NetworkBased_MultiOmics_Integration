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
library(data.table)


## 0. Definition of input/output paths ----------------------------------------------------------------------
# Input/Output path 
initial_path <- file.path("WGCNA_mRNA/data")
out_path <- file.path("WGCNA_mRNA/output")
plot_path <- file.path("WGCNA_mRNA/plot")

# # load output
bwnet <- readRDS(file.path(out_path, "bwnet.rds"))
str(bwnet)

# 2. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs
# Print out a preview
head(module_eigengenes)

# get number of genes for each module
num.per.module <- table(bwnet$colors)
    #   black        blue       brown       green greenyellow        grey 
    #     658        2945        1625        1216          33         562 
    # magenta        pink      purple         red         tan   turquoise 
    #     199         265          74         854          29        3758 
    #  yellow 
    #    1351  

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
rep('VPA_SIM', 5), rep('VS_GEMTAX', 6)))
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
nSample_tot <- nrow(norm.counts)
nGenes_tot <- ncol(norm.counts)

# qui lo calcoliamo ma in realtà la funzione della heatmap calcola la correlazione e i pvalue direttamente
# compute correlation between eigengenes of modules and traits (binarized treatments)
module.trait_tot.corr <- cor(module_eigengenes, traits_tot, use = 'p')
module.trait_tot.corr.pvals <- corPvalueStudent(module.trait_tot.corr, nSample_tot)

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
norm.counts <- readRDS(file.path(out_path, "norm_counts.rds"))
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
# (per selezionare tutti quelli con almeno 2 sterischi perchè quelli con 1 
# è il confronto che ha già altri moduli molto più significativi)
pval_star_threshold <- 0.01
significant_rows <- rownames(module.trait.corr.pvals)[apply(module.trait.corr.pvals < pval_star_threshold, 1, any)]
# Rimuove il prefisso "ME" dai nomi
colors <- sub("^ME", "", significant_rows)
colors
saveRDS(colors, file = file.path(out_path, "colors.rds"))

# Per ogni riga significativa, trova le colonne con p-value < pval_star_threshold
significant_cols <- lapply(significant_rows, function(row) {
  which(module.trait.corr.pvals[row, ] < pval_star_threshold)
})
names(significant_cols) <- colors
significant_cols

# Load gene_name, gene_id to save the lists with the gene symbol and not the ensembl
# Load txt file created by nf-core pipeline, with the number of reads for each peak 
file = file.path(initial_path, 'salmon.merged.gene_counts.tsv')
# Read the file and skip the first line + leggi solo le colonne "gene" e "count"
counts_data <- fread(file, sep = "\t", select = c("gene_id", "gene_name"))
head(counts_data)
# NB Since we want to use the gene_name column but we have repetitions, we need to change repeted names with the corresponding ensemble ids
# it is necesserary because the modules were created using the gene_id column, we can keep the name only if it is in a 1:1 relation with the id
duplicati <- names(which(table(counts_data$gene_name) > 1))

# # replace duplicated names with the content of gene_id
# for (gene in duplicati) {
#   idx <- which(counts_data$gene_name == gene)
#   counts_data$gene_name[idx] <- counts_data$gene_id[idx]
# }

# merge of module.gene.mapping to obtain the gene symbols
head(module.gene.mapping)
module.gene.mapping$gene_id <- rownames(module.gene.mapping)
merged_data <- merge(module.gene.mapping, counts_data, by = "gene_id", all.x = TRUE)
dim(module.gene.mapping)[1] == dim(merged_data)[1]
# Controlla se gene_name ha NA
print(sum(is.na(merged_data$gene_name)))
# create combined name
merged_data$combined_name <- paste0(merged_data$gene_name, "_", merged_data$gene_id)


# Creazione e salvataggio delle liste di geni in un unico ciclo
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
  # Salvataggio immediato
  writeLines(unique(gene_lists[[color]]), file.path(out_path, paste0(color, "_genes.txt")))
}

# print dimentions of modules
lapply(gene_lists, length)
        # $tan
        # [1] 29

        # $pink
        # [1] 265

        # $red
        # [1] 854

        # $greenyellow
        # [1] 33

        # $green
        # [1] 1216

        # $magenta
        # [1] 199

# 3B. Intramodular analysis: Identifying driver genes ---------------

# merge con gene_names to have gene symbols in the plots later
# Crea un nome combinato
counts_data$combined_name <- paste0(counts_data$gene_name, "_", counts_data$gene_id)
t_norm.counts <- t(norm.counts)
# da matrice a data frame
t_norm_df <- as.data.frame(t_norm.counts)
t_norm_df$gene_id <- rownames(t_norm.counts)
# merge per ottenere combined_name
merged_df <- merge(t_norm_df, counts_data[, c("gene_id", "combined_name")], 
                   by = "gene_id", all.x = TRUE)
# usare combined_name come nuovi rownames
rownames(merged_df) <- merged_df$combined_name
# rimuovi le colonne extra
merged_df$gene_id <- NULL
merged_df$combined_name <- NULL
# trasposta
norm.counts <- t(as.matrix(merged_df))

# compute correlation between eigengenes of modules and normalized counts
# Supponiamo di aver fatto il merge tra gene_id e gene_name
# merged_data è il risultato che ha sia gene_id che gene_name

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

gene_groups_comb_names <- lapply(names(significant_cols), function(color) {
  list(genes = gene_lists_comb_names[[color]], conditions = significant_cols[[color]])
})
# Assegna i nomi ai gruppi di geni
names(gene_groups_comb_names) <- names(significant_cols)


### Create heatmap for each ME module ###
#(qui consideriamo i geni significativi per almeno uno dei tratti selezionati)
for (group_name in names(gene_groups_comb_names)) {
  
  # Prendi i geni e le condizioni
  genes <- gene_groups_comb_names[[group_name]]$genes
  conditions <- gene_groups_comb_names[[group_name]]$conditions
  
  # Calcola la matrice di correlazione
  gene.signf.corr <- cor(norm.counts[, genes], traits, use = 'p')
  
  # Calcola la matrice di p-value con corPvalueStudent
  gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
  
  # Correggi i p-value per colonna mantenendo la struttura originale
  gene.signf.corr.padj <- apply(gene.signf.corr.pvals, 2, function(x) p.adjust(x, method = "BH"))

  p_threshold <- 0.01
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


### Create txt files with gene names for each ME module and each significant condition ###
#(if a module is sign for two different condition, we will have 2 txt) 
for (group_name in names(gene_groups_comb_names)) {
  
  # Prendi i geni e le condizioni
  genes <- gene_groups_comb_names[[group_name]]$genes
  conditions <- gene_groups_comb_names[[group_name]]$conditions
  
  # Calcola la matrice di correlazione
  gene.signf.corr <- cor(norm.counts[, genes], traits, use = 'p')
  # Calcola la matrice di p-value con corPvalueStudent
  gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
  # Correggi i p-value per colonna mantenendo la struttura originale
  gene.signf.corr.padj <- apply(gene.signf.corr.pvals, 2, function(x) p.adjust(x, method = "BH"))
  p_threshold <- 0.01
  
  # Seleziona i geni con padj < p_threshold nelle condizioni specificate
  for (i in seq_along(conditions)) {
    condition <- conditions[[i]]  # Ottieni la condizione
    name <- names(conditions)[i]  # Ottieni il nome della condizione
    
    selected_genes <- rownames(gene.signf.corr.padj)[which(gene.signf.corr.padj[, condition] < p_threshold)]
    # Rimuove tutto prima del primo underscore, tenendo solo il gene id (ENS)
    selected_genes_clean <- sub(".*_", "", selected_genes)
    
    # Salva la lista in un file esterno solo se contiene geni
    if (length(selected_genes_clean) > 0) {
      cat(paste0('\n', group_name, ' - ', name, ': ', length(selected_genes_clean)))
      file_name <- file.path(out_path, paste0("significant_genes_", group_name, "_", name, ".txt"))
      writeLines(selected_genes_clean, file_name)
    }
  }
}
