####################################################
### Over Representation Analysis of ours modules ###
####################################################

#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#install.packages("wordcloud")
library(clusterProfiler)
# library(wordcloud)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)
#install.packages("devtools")
#install.packages("cli")
#devtools::install_github("bzhanglab/sumer")
library(sumer)

## 1. Input and output path from $WORK
initial_path <- file.path("ShortPaper_CIBB25/data")
in_path_ora <- file.path("ShortPaper_CIBB25/WGCNA_mRNA/ORA/GSEA_MSigDB_gmt")
out_path <- file.path("ShortPaper_CIBB25/WGCNA_mRNA/ORA/output")

# Load data
# Counts data in order to obtain gene symbols
file = file.path(initial_path, 'salmon.merged.gene_counts.tsv')
counts_data <- read.delim(file, sep = "\t", header = TRUE)
dim(counts_data) #62710    36
# Filtering based on counts
smallestGroupSize <- 6 
keep <- rowSums(counts_data >= 10) >= smallestGroupSize
counts_data_filter <- counts_data[keep,]
dim(counts_data_filter) #19303
# Get gene symbols
colnames(counts_data_filter)
gene_symbol <- counts_data_filter$gene_name
length(gene_symbol) #19303
# remove ENSBL genes
gene_symbol_noENSG <- gene_symbol[!grepl("^ENSG", gene_symbol)]
length(gene_symbol_noENSG) #16243


############# ORA analysis for all the VSGT modules in a for cycle #############
in_path_genes <- 'ShortPaper_CIBB25/WGCNA_mRNA/output'
file_list <- c(
  'significant_genes_green_data.VS_GEMTAX.vs.all.txt',
  'significant_genes_magenta_data.VS_GEMTAX.vs.all.txt',
  'significant_genes_pink_data.VS_GEMTAX.vs.all.txt',
  'significant_genes_red_data.VS_GEMTAX.vs.all.txt'
)

gmt_file <- read.table(
  file.path(in_path_ora, 'GOKeggReactomeWP_sets.gmt'),
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)

############### Resolve pathway names that are too long ############
# Modify pathway name that is too long in gmt file
file_path <- file.path(in_path_ora, 'GOKeggReactomeWP_sets_INLINE.gmt')
lines <- readLines(file_path)
# Define the pattern
pattern <- "reactome_factors_involved_in_megakaryocyte_development_and_platelet_production"
new_pattern <- 'r_factors_involved_in_megakaryocyte_develop_and_platelet_prod'
# Find the line(s) containing the pattern
idx <- grep(pattern, lines)
# Replace the old pattern with the new one in those lines
lines[idx] <- sub(pattern, new_pattern, lines[idx])
# write back to the file
writeLines(lines, file_path)

for (file in file_list) {
  # Extract module color
  module_color <- sub("significant_genes_(.*?)_data.*", "\\1", file)
  cat(module_color)
  
  # Load gene list
  module_genes <- readLines(file.path(in_path_genes, file))
  cat(length(module_genes))
  
  # Ensembl to SYMBOL
  mapped_genes <- mapIds(
    org.Hs.eg.db,
    keys = module_genes,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  mapped_genes <- na.omit(mapped_genes)
  cat(length(mapped_genes))
  
  # Enrichment analysis
  enrich_result <- enricher(
    gene = mapped_genes,
    TERM2GENE = gmt_file,
    universe = gene_symbol_noENSG,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500
  )
  
  enrich_df <- data.frame(enrich_result)
  
  enrich_df$ID <- tolower(enrich_df$ID)
  
  # Modify pathway name that is too long before writing txt file
  # Replace pattern with new pattern in the ID column
  enrich_df$ID <- sub(pattern, new_pattern, enrich_df$ID)
  
  if (nrow(enrich_df) > 0) {
    enrich_df$qvalue <- -log10(enrich_df$qvalue)
    
    # Save table
    out_txt <- file.path(out_path, paste0('rna_', module_color, '_score.txt'))
    write.table(
      enrich_df[c('ID', 'qvalue')],
      file = out_txt,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )

  } else {
    message(paste("No enrichment results for module:", module_color))
  }
}


############# sumer #############
sumer("ShortPaper_CIBB25/WGCNA_mRNA/ORA/sumer_config.json", "ShortPaper_CIBB25/WGCNA_mRNA/ORA/output/sumer")
