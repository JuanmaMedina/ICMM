## DESeq2 Pipeline ##
# Why so many genes DE expression https://support.bioconductor.org/p/53957/

setwd("RNAseq/data_analysis_folder/count_data/")

# Library
library("DESeq2")
library("ggplot2")
library("dplyr")
library("pheatmap")
library("RColorBrewer") # For the hierarchical clustering coloring
library("BiocParallel")

# Reading the read-counts from each sample (FeatureCounts output)

# Controls
control_132 <- read.table("control_H132_counts_simple.txt", h = T)
control_268 <- read.table("control_H268_counts_simple.txt", h = T)

# Translocated
trans_16952 <- read.table("control_H16952_counts_simple.txt", h = T)
trans_17002 <- read.table("control_H17002_counts_simple.txt", h = T)

# Re-name the samples
names(control_132) <- c("GENE", "COUNT")
names(control_268) <- c("GENE", "COUNT")
names(trans_16952) <- c("GENE", "COUNT")
names(trans_17002) <- c("GENE", "COUNT")

## Design of the COUNT DF ##
expression_matrix <- cbind(control_132, control_268$COUNT, trans_16952$COUNT, trans_17002$COUNT)

# Genes as rownames
cts <- expression_matrix[, -1]
rownames(cts) <- expression_matrix[, 1]

# Correcting colnames
colnames(cts) <- c("control_132", "control_268", "trans_16952", "trans_17002")

# Metadata
samples <- colnames(cts)
condition <- c("control", "control", "translocation", "translocation")
type <- c("W-42", "M-57", "M-58", "W-51")

## Design of the METADATA DF ##
coldata <- data.frame(samples, condition, type)

# Samples as rownames
rownames(coldata) <- coldata[, 1]
coldata <- coldata[, -1]

# Sanity check: columns of cts and rows of coldata have to be the SAME and be in the SAME ORDER
all(colnames(cts) %in% rownames(coldata))
all(colnames(cts) == rownames(coldata))

# Construction of trh DESeqDataSet object:
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# Pre-filtering: remove rows in which there are no reads
dds <- dds[rowSums(counts(dds)) >= 1, ]

# Set factor levels (which level represents the control group and the translocated one)
dds$condition <- factor(dds$condition, levels = c("control", "translocation"))

## DE ANALYSIS ##
dds <- DESeq(dds)

# --> Parallelization step incorporation
register(MulticoreParam(4))

## RESULTS (DEFAULT ALPHA = 0.1!!!!, here is changed to 0.05) ##
res <- results(dds, alpha = 0.05)
res

# Information of each column of the results
mcols(res)$description

# Summarization of the results
summary(res)

# Number of significantly expressed genes (adjusted p-values less than 0.05) ADJUST P VALUE ACCORDING TO RES!!!
sum(res$padj < 0.05, na.rm = TRUE)

# Number of total genes
nrow(dds)


## DATA EXPLORATION ##

# MA PLOT
# Log2 FC attributables to a given variable over the mean of normalized counts for all the samples
# (red dots: adjusted p-values < 0.05; triangled points: points falling out of the window)
plotMA(res, ylim = c(-2, 2))

# Customized implementation of the MA plot
plot_res <- res
plot_res$significant = as.factor((plot_res$padj < .05))
plot_res$significant[is.na(plot_res$significant)] = F

ggplot(as.data.frame(plot_res), aes(x = log2(baseMean), y = log2FoldChange, color = significant)) +
  geom_point() +
  geom_hline(color = "blue3", yintercept = 0) +
  stat_smooth(se = FALSE, method = "loess", color = "red3") +
  scale_color_manual(values = c("Black", "Red"))


# MA PLOT (of shrunken log2 FC, obtained with a bayesian prior assumption that removes the noise
# associated with log2 FC from low count genes without requiring arbitrary filtering thresholds)
resLFC <- lfcShrink(dds, coef = 2, res = res)
plotMA(resLFC, ylim = c(-2, 2))


# Customized implementation of the shrunken MA plot
plot_res_LFC <- resLFC
plot_res_LFC$significant = as.factor((plot_res_LFC$padj < .05))
plot_res_LFC$significant[is.na(plot_res_LFC$significant)] = F

ggplot(as.data.frame(plot_res_LFC), aes(x = log2(baseMean), y = log2FoldChange, color = significant)) +
  geom_point() +
  geom_hline(color = "blue3", yintercept = 0) +
  stat_smooth(se = FALSE, method = "loess", color = "red3") +
  scale_color_manual(values=c("Black","Red"))


# PLOT COUNTS
# Count of reads for one single gene across the groups, normalizing counts by sequencing depth and adding a 
# pseudocount of 1/2 to allow for log scale plotting. In this case, the gene with smallest p-value
plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")

d <- plotCounts(dds, gene = which.min(res$padj), intgroup = "condition", returnData = TRUE)

ggplot(d, aes(x = condition, y = count)) + 
  geom_point(position = position_jitter(w = 0.1, h = 0)) + 
  scale_y_log10(breaks = c(25, 100, 400, 1600))


# Regularized log transformation to remove the dependence of the variance on the mean, particularly the high variance 
# of the logarithm of count data when the mean is low
rld <- rlog(dds, blind = FALSE)

# Variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)


# HIERARCHICAL CLUSTERING: overview over similarities and dissimilarities between samples

# Create distance matrix
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep = "-")
colnames(sampleDistMatrix) <- NULL

# Create colors (using the default ones)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# Default PCA of the regularized log transform
plotPCA(rld, intgroup = c("condition", "type"))

# Customized PCA
pcaData <- plotPCA(rld, intgroup = c("condition", "type"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition, shape = type)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed()


### -------------------------------------------------------------------- ###

# Get normalized counts
norm_dds <- estimateSizeFactors(dds)
norm_counts <- counts(norm_dds, normalized = TRUE)

# Gene of interest
rown <- c("ARHGAP42")
plotCounts(dds, gene = rown, intgroup = "condition")


log2fc_p_05 <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 2),]
log2fc_p_05 <- log2fc_p_05[order(log2fc_p_05$pvalue),]

final_res <- res[order(res$pvalue),]

# Split into up and downregulated genes
log2fc_p_05_up <- res[which(res$padj < 0.05 & res$log2FoldChange > 2), ]
log2fc_p_05_up <- log2fc_p_05_up[order(log2fc_p_05_up$pvalue), ]

log2fc_p_05_down <- res[which(res$padj < 0.05 & res$log2FoldChange < -2), ]
log2fc_p_05_down <- log2fc_p_05_down[order(log2fc_p_05_down$pvalue), ]


# Write output
# write.csv(as.data.frame(res), file="all_genes_results.csv",quote = F)






