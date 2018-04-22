# Set cd
setwd("karyotypes_niels/all_sheets_karyos/")

library("stringr")

### Regex design ### 
# This regex captures mut + chr + cytobands (e.g. t(8;19)(p10;q10))
regex_cyto <- regex("[itd].+\\)")

# This regex captures un-separated cytobands (e.g. q11.2q26.1)
regex_semico <- regex("([0-9])([a-z])")  

# This regex captures random erroneous characters (e.g. inv(Y)(p11.2p11.2?3)) 
regex_characters <- regex("[\\?~]")

# # This regex captures contiguous "." and ";" originated from the removal of erroneous chars
# regex_contig <- regex("\\.;")

### DB file processing ###

# Load UCSC database
db <- read.table("merged_db.bed")
db <- db[, c(2, 3, 4, 1)]

# Final format -> chrX:cytobandY
db$final_cyto <- apply(db[, c("V2", "V1")], 1, paste, collapse = ":")

cytoband_matcher <- function (d)  {
  
  # Load file, remove duplicates and save as total records
  dfX <- read.table(d, sep = "\n")
  dfX <- subset(dfX, !duplicated(dfX$V1))
  
  # Clone total records and start pipeline
  df <- dfX
  
  # Capture mut + chr + cytobands in second column
  df$V2 <- str_extract(dfX$V1, regex_cyto)
  
  # Remove NAs
  df <- df[complete.cases(df), ]
  
  # Keep rows with correct cytobands, formed by 4 parenthesis: 2 "(" + 2 ")" 
  df <- df[which(str_count(df$V2, "\\(") == 2), ]
  df <- df[which(str_count(df$V2, "\\)") == 2), ]
  
  # Replace ":" or "-" separated cytobands information by ";"
  df$V2 <- gsub("[:-]", ';', df$V2)
  
  # Remove blanks
  df$V2 <- gsub(" ", "", df$V2)
  
  # Separate un-separated cytobands
  df$V2 <- gsub(regex_semico, '\\1\\;\\2', df$V2)
  
  # Remove erroneous characters
  df$V2 <- gsub(regex_characters, "", df$V2)
  
  # Remove rows with no ";" or more than 2 ";" (cytoband info messed up)
  df <- df[which(str_count(df$V2, "\\;") != 0), ]
  df <- df[which(str_count(df$V2, "\\;") <= 2), ]
  
  # Replace comma separated-sub-band information by dots
  df$V2 <- gsub(',', '.', df$V2)
  
  # # Keep rows with cytoband formed by > 10 and < 24 chars
  # df <- subset(df, nchar(V2) > 10 & nchar(V2) < 24)
  
  # Rename headers
  colnames(df) <- c("-", "info")
  
  # Remove duplicated cyotbands
  df <- subset(df, !duplicated(df$info))
  
  # Split mut + chr from cytoband
  # Strsplit: splitting by ESCAPED symbol (")") over CHARACTER vector
  df$nature <- lapply(strsplit(as.character(df$info), "\\)"), "[", 1)
  df$cytoband <- lapply(strsplit(as.character(df$info), "\\)"), "[", 2)
  
  # Remove initial "("
  df$cytoband <- substring(df$cytoband, 2)
  
  # Split the two cytobands involved in the mutation
  df$cyto_1 <- lapply(strsplit(as.character(df$cytoband), "\\;"), "[", 1)
  df$cyto_2 <- lapply(strsplit(as.character(df$cytoband), "\\;"), "[", 2)
  
  # Remove patients with no cytoband information
  df <- na.omit(df[, -1])
  
  # Split chromosomes from mutation type (inversions affect 1 chr, translocations affect 2 chr)
  df$type <- lapply(strsplit(as.character(df$nature), "\\("), "[", 1)
  df$chr <- lapply(strsplit(as.character(df$nature), "\\("), "[", 2)
  
  # Split chrA from chrB
  df$chrA <- lapply(strsplit(as.character(df$chr), "\\;"), "[", 1)
  df$chrB <- lapply(strsplit(as.character(df$chr), "\\;"), "[", 2)
  
  # For inversions, with NA in "chrB", duplicate the chromosome appearing in "chrA"
  df$chrB[is.na(df$chrB)] <- as.character(df$chrA[is.na(df$chrB)])
  
  # Add "chr" to "chrA" and "chrB" columns
  df$chrA <- paste("chr", df$chrA, sep = "")
  df$chrB <- paste("chr", df$chrB, sep = "")
  
  # Paste chr to cytoband (final format -> chrX:cytoband)
  df$final_cyto_A <- apply(df[, c("chrA", "cyto_1")], 1, paste, collapse = ":")
  df$final_cyto_B <- apply(df[, c("chrB", "cyto_2")], 1, paste, collapse = ":")
  
  # Final DFs
  df <- df[, c("info", "final_cyto_A", "final_cyto_B")]
  
  # Remove sub-band information
  df$final_cyto_A <- sub("\\..+", "", df$final_cyto_A)
  df$final_cyto_B <- sub("\\..+", "", df$final_cyto_B)
  
  # Matching with DB
  df$coord1 <- db[match(df$final_cyto_A, db$final_cyto), ]
  df$coord2 <- db[match(df$final_cyto_B, db$final_cyto), ]
  
  # # Record no-matching breakpoints
  # df2 <- df[rowSums(is.na(df[ , c(4, 5)])) > 0, c(1:3)]
  
  # Remove rows without two matches in the DB (i.e. non identified BP)
  df <- df[rowSums(is.na(df[ , c(4, 5)])) == 0, ]
  
  # Record no matching breakpoints
  df2 <- as.data.frame(dfX[rownames(dfX)[which(!rownames(dfX) %in% rownames(df))], ])
  
  # Select columns
  df <- cbind(df$info, 
              df[, c("coord1")][c("V2", "V3", "V4", "final_cyto")], 
              df[, c("coord2")][c("V2", "V3", "V4", "final_cyto")])
  
  # Rename
  colnames(df) <- c("info", 
                    "chrA", "stA", "endA", "infoA", 
                    "chrB", "stB", "endB", "infoB")
  
  # Save DF: only chr, start & end
  df[, c(2,3,4,6,7,8)]
  
}

# carrier_familial_1 <- cytoband_matcher("carrier_familial_1.txt")
# abortions_2 <- cytoband_matcher("abortions_2.txt")
# affected_3 <- cytoband_matcher("affected_3.txt")
# prenatal_carriers_4 <- cytoband_matcher("prenatal_carriers_4.txt")
# prenatal_affected_5 <- cytoband_matcher("prenatal_affected_5.txt")
# indication_not_pre_6 <- cytoband_matcher("indication_not_pre_6.txt")
# indication_not_post_7 <- cytoband_matcher("indication_not_post_7.txt")
# post_mosaic_8 <- cytoband_matcher("post_mosaic_8.txt")

controls <- cytoband_matcher("juanma_controls.txt")
affected <- cytoband_matcher("juanma_affected.txt")

# # Writing
# write.table(x = df2[, 1],
#             file = "C:/Users/Juanma/Desktop/karyotypes_niels/nomatch_controls.txt",
#             quote = F, row.names = F, col.names = F)

# My total number of hits

# ff <- list(abortions_2, affected_3, carrier_familial_1, indication_not_post_7,
#            indication_not_pre_6, post_mosaic_8, prenatal_affected_5, prenatal_carriers_4)
# 
# sum(sapply(ff, nrow))


# RCircos plots

library("RCircos")

## Load chromosome ideogram data: 3 columns ##
# Chr. names, start and end
data("UCSC.HG19.Human.CytoBandIdeogram")
head(UCSC.HG19.Human.CytoBandIdeogram)

## Breakpoints ##

# Each rows is a pair of breakpoints: 6 columns
# Chr. A names, chr. A start and chr. A end; chr. B names, chr. B start and chr. B end
# my_links <- rbind(carrier_familial_1, abortions_2, affected_3, prenatal_carriers_4, 
#                   prenatal_affected_5, indication_not_pre_6, indication_not_post_7, post_mosaic_8)

my_links <- rbind(affected, controls)

# REMOVE DUPLICATES HERE!

## TADs represented as vertical lines: 3 columns ##
# Chr. names, start and end
tad <- read.table(file = "../../../Mana/Mana_tads_bp/Circ_plot/tad_coordinates.txt")

# Initialize chr. ideogram
RCircos.Set.Core.Components(cyto.info = UCSC.HG19.Human.CytoBandIdeogram, 
                            chr.exclude = NULL, tracks.inside = 10, tracks.outside = 0)

# Building the plot

# 1. Initialize graphic device
out.file <- "plot_both.pdf"
pdf(file = out.file)
RCircos.Set.Plot.Area()

# 2. Plot chromosome ideogram
RCircos.Chromosome.Ideogram.Plot(tick.interval = 50)
RCircos.Draw.Chromosome.Ideogram()
RCircos.Label.Chromosome.Names(chr.name.pos=NULL)

## Modifiy some parameters

RCircos.List.Plot.Parameters()

rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$line.color <- "red"
RCircos.Reset.Plot.Parameters(rcircos.params)


# Links
RCircos.Link.Plot(my_links, track.num = 2, TRUE)

# Vertical lines
RCircos.Vertical.Line.Plot(tad, track.num = 1, side = "in")

dev.off()
