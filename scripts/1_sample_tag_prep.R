suppressPackageStartupMessages({
  library(ShortRead)
  library(Matrix)
  library(tidyverse)
})

ST_fasta <- readDNAStringSet(snakemake@input[["fasta"]])
ST_fasta <- subseq(ST_fasta, start = 26, end = 56)
st_names <- names(ST_fasta)
st_names <- gsub("\\|.*", "", st_names)
st_seq <- paste(ST_fasta)
df <- data.frame(st_seq, st_names)
write_tsv(df, snakemake@output[["tsv"]], col_names = FALSE)