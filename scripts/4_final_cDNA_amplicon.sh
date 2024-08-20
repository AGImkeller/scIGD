#!/bin/bash

alleles_dir="$1"
gene_names="$2"
cDNA_RNA="$3"
output_cDNA="$4"
output_lookup_table="$5"

conc_alleles=$(awk '
    /^>/ {
        if (header != "") {
            headers[header] = headers[header] typed_alleles
        }
        header = $0
        typed_alleles = ""
        next
    }
    /^Typed alleles:/ {
        getline
        while (!/^----------------------$/) {
            if (/^[[:space:]]*[0-9]+/) {
                sub(/^[[:space:]]*[0-9]+[[:space:]]*/, "")
            }
            typed_alleles = typed_alleles $0 "\n"
            getline
        }
    }
    END {
        if (header != "") {
            headers[header] = headers[header] typed_alleles
        }
        for (header in headers) {
            printf("%s\n%s", header, headers[header])
        }
    }
' "$alleles_dir"/alleles* | awk '!seen[$0]++')

cDNA_HLA=$(echo "$conc_alleles" | awk '/^>/{if (header != "") sequences=0; header=$0; sub(/\|.*/, "", header); next} {sequences++; print header "-AV-" sequences; print}')

cDNA_RNA_filtered=$(awk '(NR==FNR) { toRemove[$1]; next }
    /^>/ { p=1; for(h in toRemove) if ( h ~ $0) p=0 }
    p' "$gene_names" "$cDNA_RNA")

echo -e -n "${cDNA_RNA_filtered}\n${cDNA_HLA}" > "$output_cDNA"

echo "Allele,Gene,Function" > "$output_lookup_table" && echo "$cDNA_HLA" | awk '/^>/ { header=substr($0,2); gene=substr(header,1,index(header,"-AV-")-1); print header "," gene }' | awk -F, 'BEGIN {OFS=","} NR>0 {if ($2=="HLA-A" || $2=="HLA-B" || $2=="HLA-C") $3="HLA_class_I"; else if ($2=="HLA-DRA" || $2=="HLA-DRB1" || $2=="HLA-DRB3" || $2=="HLA-DRB5" || $2=="HLA-DQA1" || $2=="HLA-DQB1" || $2=="HLA-DPA1" || $2=="HLA-DPB1") $3="HLA_class_II"; else if ($2=="HLA-E" || $2=="HLA-F" || $2=="HLA-G" || $2=="HLA-DMA" || $2=="HLA-DMB" || $2=="HLA-DOA" || $2=="HLA-DOB") $3="HLA_non_classical"; else if ($2=="HLA-H" || $2=="HLA-J" || $2=="HLA-K" || $2=="HLA-L") $3="HLA_non_coding"} 1' >> "$output_lookup_table"