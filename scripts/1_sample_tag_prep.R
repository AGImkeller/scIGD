# R

## this R script reads the Sample Tags FASTA file and generates a TSV file containing the sequences (1st column) and their names (2nd column)
## this TSV file is necessary for running kallisto with kite workflow for sample tag quantification

suppressPackageStartupMessages({
  library(ShortRead)
  library(Matrix)
  library(tidyverse)
})

ST_fasta <- readDNAStringSet(snakemake@input[["sample_tag_fasta"]])

## subsetting the sample tag sequences to 31bp
## kallisto's upper limit for its pseudo-alignment process with kite workflow is 31bp
## kite workflow is designed for feature barcoding -- recommended to use for sample tag quantification
## since demultiplexing is implemented specifically for BD Rhapsody data
## the chosen start and end bases allow for proper differentiation between BD's sample tag sequences
ST_fasta <- subseq(ST_fasta, start = 26, end = 56)

st_names <- names(ST_fasta)

## shortening the sample tag names 
## this is done simply to make them easier to read by eye
st_names <- gsub("\\|.*", "", st_names)

st_seq <- paste(ST_fasta)
df <- data.frame(st_seq, st_names)
write_tsv(df, snakemake@output[["sample_tag_tsv"]], col_names = FALSE)