suppressPackageStartupMessages({
  library(ShortRead)
  library(Matrix)
  library(tidyverse)
  library(sparseMatrixStats)
})

dir.tags <- snakemake@input[["dir"]]
tags <- readMM(paste0(dir.tags, "/counts_unfiltered/cells_x_features.mtx", ""))
cells <- read.table(paste0(dir.tags, "/counts_unfiltered/cells_x_features.barcodes.txt", ""), header = FALSE)
samples <- read.table(paste0(dir.tags, "/counts_unfiltered/cells_x_features.genes.txt", ""), header = FALSE)
rownames(tags) <- cells$V1
colnames(tags) <- samples$V1
tags <- as(tags, "CsparseMatrix")
df1 <- data.frame(Cellular.Barcode = rownames(tags),
                  Sample.Tag = ifelse(sparseMatrixStats::rowMaxs(tags) >= 0.75 * rowSums(tags),
                                      colnames(tags)[max.col(tags)],
                                      "Multiplet"))
df1 <- df1[!df1$Sample.Tag == 'Multiplet', ]
write.table(df1, file = snakemake@output[["df1"]], 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)