## DESEq2-GAGE Pipeline ##

setwd("/media/bioinfor1/9c6bddae-3750-4442-9871-6bfcd83cc0b2/home/binf1/inigo/manna_rna_seq/Mana_RNAseq/data_analysis_folder/count_data/")

# Library
library("DESeq2")
library("gage")
library("gageData")                            # Pre-compiled DB mapping genes to KEGG pathways and GO terms for common organisms
library("pathview")                            # Draw pathway diagrams
library("biomaRt")


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

## RESULTS (change default alpha to 0.05) ##
res <- results(dds, alpha = 0.05)

## GAGE analysis ##

## KEGG pathway analysis ##
# DE genes as probe set: APPLY A SECOND CRITERIA OF GENES LFC WITH ADJUSTED PVAL < 0.05 (from 28872 to 1965)
my_res <- data.frame(res[c(2, 6)])
my_res <- subset(my_res, padj <= 0.05, select = log2FoldChange) 

# Rownames as first column, and change its name to merge in a posterior step
my_res <- data.frame(names = row.names(my_res), my_res)
names(my_res)[1] <- "SYMBOL"

# Entrez IDs for my probe, with some NAs corresponding to not found genes
id.map.symbol <- id2eg(ids = rownames(my_res), category = "SYMBOL", org = "human")

# Merge Entrez IDs and LFC by symbol, remove NAs and reorder 
ids <- merge(my_res, id.map.symbol, by = "SYMBOL")
ids <- data.frame(ids[complete.cases(ids), ])

# Keep Entrez ID and LFC, and set ID as rownames
ids <- ids[, c(3, 2)]
my_ids <- data.frame(ids[, -1], row.names = ids[, 1])

# GAGE analysis for DE genes
data("kegg.gs")
fc.kegg.p <- gage(exprs = my_ids, gsets = kegg.gs, ref = NULL, samp = NULL)

# GAGE pathway enrichment analysis consider only LFC and not the number of samples that went into the LFC comparison, which
# causes q-values from GAGE to be falsely inflated. POTENTIAL SOLUTION: DE table filtered with a p-val < 0.05 AND FDR < 0.01,
# so only significantly expressed genes with a low FDR are used in the calculation of enriched pathways. All KEGG pathways IDs
# with a p-value < 0.05 are returned

# Article: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4668955/

# But in this example, only upreg genes have inflated qvals, the downreg do not!!!! <----------------

upreg_genes <- data.frame(fc.kegg.p$greater)
downreg_genes <- data.frame(fc.kegg.p$less)

# List of up and down-regulated pathways
head(upreg_genes[order(upreg_genes$q.val), ], 10)
head(downreg_genes[order(downreg_genes$q.val), ], 10) 

# List of up and down-regulated pathways selected based on q-val < 0.1
upreg_genes[upreg_genes$q.val < 0.1 & !is.na(upreg_genes$q.val), ]
downreg_genes[downreg_genes$q.val < 0.1 & !is.na(downreg_genes$q.val), ]

# List of up-regulated pathways selected based on q-val < 0.1
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 &
  !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]

# List of down-regulated pathways selected based on q-val < 0.1
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 &
  !is.na(fc.kegg.p$less[, "q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]

# IDs of all the pathways, UPr and DOWNr
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)

# DOWNLOAD the pathways extracted
pv.out.list <- sapply(path.ids2[1:2], function(pid) 
  pathview(gene.data =  my_ids, pathway.id = pid,
           species = "hsa", out.suffix = "deseq2"))


## GO analysis ##
# Load list of 17000 elements, each of this is a character vector of Entrez IDs for a single GO term
data("go.sets.hs") 

# Load named list of three elements: biological processes, cellular components and molecular functions. I will run a GO
# enrichment for each of the three subtrees
data("go.subs.hs")

my_ids.bp.p <- gage(exprs = my_ids, gsets = go.sets.hs[go.subs.hs$BP], ref = NULL, samp = NULL)
my_ids.mf.p <- gage(exprs = my_ids, gsets = go.sets.hs[go.subs.hs$MF], ref = NULL, samp = NULL)
my_ids.cc.p <- gage(exprs = my_ids, gsets = go.sets.hs[go.subs.hs$CC], ref = NULL, samp = NULL)

# List of Biological Processes GO terms for up and down-regulated genes
upreg_genes_bp <- data.frame(my_ids.bp.p$greater)
downreg_genes_bp <- data.frame(my_ids.bp.p$less)

head(upreg_genes_bp[order(upreg_genes_bp$q.val), ], 10)
head(downreg_genes_bp[order(downreg_genes$q.val), ], 10) 

# List of Molecular Functions GO terms for up and down-regulated genes
upreg_genes_mf <- data.frame(my_ids.mf.p$greater)
downreg_genes_mf <- data.frame(my_ids.mf.p$less)

head(upreg_genes_mf[order(upreg_genes_mf$q.val), ], 10)
head(downreg_genes_mf[order(downreg_genes_mf$q.val), ], 10) 

# List of Celular Components GO terms for up and down-regulated genes
upreg_genes_cc <- data.frame(my_ids.cc.p$greater)
downreg_genes_cc <- data.frame(my_ids.cc.p$less)

head(upreg_genes_cc[order(upreg_genes_cc$q.val), ], 10)
head(downreg_genes_cc[order(downreg_genes_cc$q.val), ], 10)

# ----------------------------------------------------------------------------------------------------- #


# Search for the GO terms associated with the interest gene: ARHGAP42
# Entrez Gene: 143872
# Ensembl: ENSG00000165895

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
mart_results <- getBM(attributes = c("refseq_dna", "go_molecular_function_id",
                                "go_molecular_function"), filters = "refseq_dna",
                 values = c("ENSG00000165895"), mart = mart)

listAttributes(
)



# Test for over-representation of gene ontology (GO) terms or KEGG pathways in a -ID transformed- Entrez Gene IDs vector
library("limma")
my_go <- goana(de = as.vector(id.map.symbol), coef = deseq2.fc, geneid = rownames(id.map.symbol), FDR = 0.1, 
               species = "Hs", species.KEGG = "hsa", gene.pathway = getGeneKEGGLinks(species.KEGG),
               pathway.names = getKEGGPathwayNames(species.KEGG, remove = TRUE, plot = TRUE))

# Results
head(my_go[order(my_go$P.DE), ], 50)


write.table(as.data.frame(upreg_genes), file="GO/UP_genes_KEGG_pathways.txt", quote = F)
write.csv(as.data.frame(downreg_genes), file="GO/DOWN_genes_KEGG_pathways.csv", quote = F)


write.

# Transformation of bioMART output to .csv
down_go_terms <- read.table("GO/down_goterms_final.txt",  header = T, sep = "\t")
