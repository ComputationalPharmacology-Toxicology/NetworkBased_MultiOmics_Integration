## FIRST PART OF WGCNA: create WGCNA object to obtain modules

# # install libraries
# BiocManager::install("impute") # in order to install WGCNA
# BiocManager::install("WGCNA") # in order to install WGCNA
# install.packages('flashClust')
# BiocManager::install("rversions") # in order to install devtools
# install.packages("xml2") # in order to install devtools
# install.packages("devtools")
# BiocManager::install('devtools')
# devtools::install_github("kevinblighe/CorLevelPlot")

# load libraries
library(devtools)
library(WGCNA)
library(flashClust)
library(curl)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(gridExtra)
library(CorLevelPlot)
library(pheatmap)


## 0. Definition of input/output paths ----------------------------------------------------------------------
initial_path <- file.path("WGCNA_miRNA/data")
out_path <- file.path("WGCNA_miRNA/output")
plot_path <- file.path("WGCNA_miRNA/plot")

## 1. Get data from nf-core piepeline ----------------------------------------------------------------------
# Load txt file with the number of reads
file = file.path(initial_path, 'mature_counts.csv')
counts_data <- read.csv(file, sep = ",", header = TRUE, row.names=1)
counts_data <- counts_data[order(rownames(counts_data)), ] # order samples by name
rownames(counts_data)
counts_data <- t(counts_data)

dim(counts_data)
head(counts_data,10)

# Check if all number are integers 
all(sapply(counts_data, is.integer))


## 2A. Detect & exclude outliers ----------------------------------------------------------------------
gsg <- goodSamplesGenes(t(counts_data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# exclude genes
data <- counts_data[gsg$goodGenes == TRUE, ]

# check again samples for outliers with another method (clustering)
htree <- hclust(dist(t(data)), method='average')
# plot to see if there are outliers
pdf(file.path(plot_path, '1_htree_out_detect.pdf'))
plot(htree)
dev.off()

# check again samples for outliers with another method (pca)
pca <- prcomp(t(data))
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
pca.dat <- as.data.frame(pca.dat)
# plot to see if there are outliers
pdf(file.path(plot_path, '1_pcaplot_out_detect.pdf'))
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
dev.off()


# 2B. Exclude genes with low variance ----------------------------------------------------------------------
# We are skipping this step for miRNAs (only did it for mRNAs)
# geneVariance <- apply(data, 1, var)
# threshold <- quantile(geneVariance, 0.10)  
data_filtered <- data#[geneVariance > threshold, ]

dim(data)
dim(data_filtered)


# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset
colnames(data_filtered)
# Create colData
sampletype <- as.factor(c(rep('CTR',6), rep('GEMTAX', 5), rep('SIM',6), rep('VPA', 6), 
rep('VPASIM', 5), rep('VSGT', 6)))
colData <- data.frame(sampletype, row.names = colnames(data_filtered))

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data_filtered))
all(rownames(colData) == colnames(data_filtered))

# create dds
dds <- DESeqDataSetFromMatrix(countData = data_filtered,
                              colData = colData,
                              design = ~ 1) # not spcifying model

dim(data_filtered)
## remove all genes with counts < 1 in more than 25% of samples (34*0.25=8.5)
## less stringent than what suggested by WGCNA on RNAseq FAQ because here we have miRNAs
dds75 <- dds[rowSums(counts(dds) >= 1) >= 8,]
nrow(dds75) #1419 genes 

# perform variance stabilization
# dds_norm <- vst(dds)
dds_norm <- varianceStabilizingTransformation(dds75) # having few rows we can not use vst()

# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()
# save for later
saveRDS(norm.counts, file = file.path(out_path, "norm_counts.rds"))


# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers to then select the best
power <- c(c(1:10), seq(from = 12, to = 30, by = 2))
# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)

sft.data <- sft$fitIndices

# visualization to pick power
pdf(file.path(plot_path, '1_rsqd_meank_threshold.pdf'))
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()
a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 100, color = 'red') +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
grid.arrange(a1, a2, nrow = 2)
dev.off()

# based on the plots we chose 6
soft_power <- 6

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

# compute adjacency
adjacency <- adjacency(norm.counts, 
                       type = "signed",
                       power = soft_power)
str(adjacency)
# save for later
saveRDS(adjacency, file = file.path(out_path, "adjacency.rds"))

# we need to set it before running blockwiseModules()
temp_cor <- cor
cor <- WGCNA::cor

# create modules
bwnet <- blockwiseModules(norm.counts,
                 maxBlockSize = nrow(dds75)+10,
                 TOMType = "signed",
                 networkType = "signed",
                 saveTOMs = FALSE,
                 power = soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,
                 verbose = 3)

cor <- temp_cor

# save output
saveRDS(bwnet, file = file.path(out_path, "bwnet.rds"))
