#!/bin/bash

fastq="$1"
seq="$2"
num_processes="$3"

ls "$fastq"/*.R2.fastq.gz | xargs -I {} -P "$num_processes" sh -c 'file="{}"
    output_file="data/allele_typing/alleles/alleles_$(basename {} .R2.fastq.gz).txt"
    while IFS= read -r header; do

        read -r sequence
        match=$(zcat "$file" | grep "^$sequence")
        sum_counts=$(echo "$match" | wc -l)
        allele_counts=$(echo "$match" | grep -v "^$" | sed -e "s/,$//" -e "s/.*\.\(.*\)\.\(.*\)$/\1.\2/" | sort | uniq -c)
        top_alleles=$(echo "$allele_counts" | sort -k1n,1 | tail -2)
        percentages=$(echo "$top_alleles" | awk -v total="$sum_counts" '\''{printf "%.2f%%\n", ($1 * 100 / total)}'\'')
        allele1_percentage=$(echo "$percentages" | awk 'NR==1{print}')
        allele2_percentage=$(echo "$percentages" | awk 'NR==2{print}')
        percentage_difference=$(awk -v a="$allele2_percentage" -v b="$allele1_percentage" '\''BEGIN {printf "%.2f\n", a - b }'\'')
        percentage_difference=$(echo "$percentage_difference" | sed 's/,/./g')
        percentages=$(echo "$percentages" | sed 's/,/./g')

        if [ "$(awk -v pd="$percentage_difference" "BEGIN {printf (pd > 50) ? "1" : "0"}")" -eq 1 ]; then
            typed_allele=$(echo "$top_alleles" | tail -1)
        else
            typed_allele=$(echo "$top_alleles")
        fi

        echo "$header: $sequence" >> "$output_file"
        echo "Top 2 alleles:" >> "$output_file"
        echo "$top_alleles" >> "$output_file"
        echo "Total number of reads: $sum_counts" >> "$output_file"
        echo "Percentage of the top 2 alleles:" >> "$output_file"
        echo "$percentages" >> "$output_file"
        echo "Difference between the top 2 alleles: $percentage_difference%" >> "$output_file"
        echo "Typed alleles:" >> "$output_file"
        echo "$typed_allele" >> "$output_file"
        echo "----------------------" >> "$output_file"

    done < data/allele_typing/short_seq.fasta
'
