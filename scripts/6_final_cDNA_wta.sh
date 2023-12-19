#!/bin/bash

allele_cDNA="$1"
genes_cDNA="$2"
gene_names="$3"
final_cDNA="$4"

transcripts=$(grep -E "gene_name:($gene_names)" "$genes_cDNA")

filtered_cDNA=$(cat "$genes_cDNA" | grep -v "$transcripts")

filtered_cDNA_1=$(echo "$filtered_cDNA" | awk '/^>/ {sub(/ gene_name.*/, ""); sub(/ gene_id:/, " "); print $0} /^[^>]/ {print}')

allele_cDNA_1=$(cat "$allele_cDNA" | sed -e 's/^>HLA:HLA\(.*\) \(.*\) \([0-9]* bp\)/>HLA\1 \2/')

cDNA="$filtered_cDNA_1
$allele_cDNA_1"

echo "$cDNA" > "$final_cDNA"