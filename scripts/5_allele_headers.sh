#!/bin/bash

## ---------------------------------------

## this script processes arcasHLA's genotype data to extract the alleles and 
## then matches them with corresponding headers in the HLA reference cDNA file 
## the goal is to create a FASTA file containing the headers of matching alleles for proper sequence extraction
## and integration with reference cDNA file later

## main input:
## - arcasHLA genotype output JSON file (referred to in this script as genotype="$1")
## - a CSV file listing all known HLA alleles (referred to in this script as list="$2")
## - a FASTA file including the DNA sequence for all known HLA alleles (referred to in this script as cdna="$3")
## - output FASTA file for matched allele headers (referred to in this script as allele_fasta="$4")

## this script extracts allele names from the genotype file and looks for matches in the specified CSV file
## then retrieves headers from the HLA reference cDNA file that correspond to the matched alleles
## to facilitate sequence extraction later

## ---------------------------------------

## input data
genotype="$1"
list="$2"
cdna="$3"
allele_fasta="$4"

## reads the arcasHLA genotype file that has a specific JSON format
## extracts the names of all alleles typed and assigns the result to 'genotype_txt'
genotype_txt=$(cat "$genotype" | grep -o '\[[^]]*\]' | sed 's/\[\|\]//g' | sed 's/ /\n/g' | tr -d ',' | tr -d '"' | uniq)

## iterates over each unique pattern from 'genotype_txt' and searches for the first matching line in "$list"
## assigns the results to 'allele_txt'
allele_txt=$(echo "$genotype_txt" | while read pattern; do cat "$list" | grep -m 1 -F "$pattern"; done)

## reads the HLA ref cDNA file, finds lines that match the allele names in 'allele_txt'
## removes the '>' character from the header and assigns the cleaned headers to 'matching_headers'
matching_headers=$(cat "$cdna" | grep -w -F "$allele_txt" | sed 's/^>//g')

echo "$matching_headers" > "$allele_fasta"