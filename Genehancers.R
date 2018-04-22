db <- read.csv(file = "genehancer.csv")
input_genes <- read.csv(file = "cneTAD_target_genes.csv", header = F)

library("stringr")
library("dplyr")
library("data.table")
library("reshape2")

# Split genhancer DB by ";"
enhancer_gene_ids_scores <- as.data.frame(str_split_fixed(db$attributes, ";", n = 160))

# Number of genes in last column specified
which(enhancer_gene_ids_scores$V160 != "")

# Keep only enhancer ID (col 1) and gene IDs (even columns)
enhancers_gene_ids <- enhancer_gene_ids_scores[, c(1, seq(2, 160, 2))]

# Remove annoying strings
ids <- data.frame(lapply(enhancers_gene_ids, function(x) {
  gsub("connected_gene=", "", x)
  }))

# Helper function to remove NAs
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

# Transform every empty cell from the DF to NAs
ids <- ids %>% mutate_each(funs(empty_as_na))

# Add chromosomic coordinates of each enhancer
ids <- data.frame(ids, db[, c(1, 4, 5)])

# Extract input genes IDs
input_genes <- as.vector(input_genes$V4)

# Enhancers that are associated at least with 1 gene of the input list
filtered_ids <- ids[which(apply(ids, 1, function(r) any(r %in% input_genes))), c(1:84)]

# Re-order columns
filtered_ids <- filtered_ids[, c(1, 82, 83, 84, 2:81)]

# Melt dataframe
melted_ids <- melt(filtered_ids, id.vars = c("V1", "chrom", "start", "end"), na.rm = T)

# Remove rows with no matching genes 
genehancers <- melted_ids[melted_ids$value %in% input_genes, c(-5)]

# Sanity checks
all(genehancers$value %in% input_genes)
sum(input_genes %in% genehancers$value)
sum(duplicated(genehancers$start) && duplicated(genehancers$end))
!any(genehancers$value %in% input_genes)

# Sorting
genehancers <- genehancers[order(genehancers$value), ]

# Remove annoying strings
genehancers <- data.frame(lapply(genehancers, function(x) {
  gsub("genehancer_id=", "", x)
}))

# Paste gene-name as suffix of each enhancer
genehancers$value <- paste(genehancers$value, genehancers$V1, sep = "_")

# Remove enhancers column
genehancers <- genehancers[, c(-1)]

# Error for hg37, it works for hg38
# genehancers <- genehancers[-(5158:5165), ]


write.table(x = genehancers, 
          file = "genes_enhancers.txt", quote = F, row.names = F, col.names = F)

