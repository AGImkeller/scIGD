#!/bin/bash

genotype="$1"
list="$2"
cdna="$3"
allele_fasta="$4"

genotype_txt=$(cat "$genotype" | grep -o '\[[^]]*\]' | sed 's/\[\|\]//g' | sed 's/ /\n/g' | tr -d ',' | tr -d '"' | uniq)

allele_txt=$(echo "$genotype_txt" | while read pattern; do cat "$list" | grep -m 1 -F "$pattern"; done)

matching_headers=$(cat "$cdna" | grep -w -F "$allele_txt" | sed 's/^>//g')

echo "$matching_headers" > "$allele_fasta"