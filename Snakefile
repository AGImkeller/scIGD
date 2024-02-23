## configuration file
configfile: 'config.yaml'

## run all analyses
rule all:
    input:
        "data/quant/kb_output_amplicon" if not config['wta'] else "data/quant/kb_output_WTA"

## ------------------------------------------------------------------------------------ ##
## demultiplexing, if needed
## ------------------------------------------------------------------------------------ ##
if config['multiplex']:
    
    ## sample tag tsv file preparation 
    rule sample_tag_prep:
        input: 
            fasta=config['sample_tag_seqs']
        output:
            tsv="data/demultiplex/meta/sample_tag.tsv"
        log:
            "logs/sample_tag_prep.log"
        script:
            "scripts/1_sample_tag_prep.R"

    ## building a sample tag index
    rule sample_tag_index:
        input:
            tsv="data/demultiplex/meta/sample_tag.tsv"
        output:
            index="data/demultiplex/meta/st_index.idx",
            fasta="data/demultiplex/meta/st_mismatch.fa",
            t2g="data/demultiplex/meta/st_f2g.txt"
        log:
            "logs/sample_tag_index.log"
        shell:
            "kb ref -i {output.index} -f1 {output.fasta} -g {output.t2g} --workflow kite {input.tsv} > {log} 2>&1"

    ## sample tag quantification
    rule sample_tag_quant:
        input:
            raw_data_fastq_list=config['raw_data_fastq_list'][0:2],
            index="data/demultiplex/meta/st_index.idx",
            t2g="data/demultiplex/meta/st_f2g.txt"
        output:
            directory("data/demultiplex/quant")
        params:
            tech=config['sequencing_technology']
        log:
            "logs/sample_tag_quant.log"
        resources:
            mem_gb=8,
            cores=4
        shell:
            "kb count -i {input.index} -g {input.t2g} -x {params.tech} -o {output} --workflow kite {input.raw_data_fastq_list} > {log} 2>&1"

    ## dataframe 1: cellular barcode + sample tag
    rule sample_tag_ident:
        input:
            dir="data/demultiplex/quant"
        output:
            df1=temp("data/demultiplex/splitting/df1.txt")
        log:
            "logs/sample_tag_ident.log"
        script:
            "scripts/2_sample_tag_ident.R"

    ## extracting cellular barcodes and their corresponding read IDs
    rule sample_tag_id:
        input:
            fastq1=config['raw_data_fastq_list'][0]
        output:
            cb=temp("data/demultiplex/splitting/cb.txt"),
            id=temp("data/demultiplex/splitting/id.txt")
        params: 
            read_len=config['read_bp_length']
        log:
            "logs/sample_tag_id.log"
        threads: 8
        run:
            shell(
                """
                zcat {input.fastq1} | awk 'NR%4==2 && length($1)=={params.read_len} {{ print substr($1, 1, 9) substr($1, 22, 9) substr($1, 44, 9) }}' > {output.cb}
                seqkit seq {input.fastq1} -m {params.read_len} -n -i > {output.id}
                """
            )

    ## dataframe 2: cellular barcode + ID
    rule sample_tag_id2:
        input:
            cb="data/demultiplex/splitting/cb.txt",
            id="data/demultiplex/splitting/id.txt"
        output:
            df2=temp("data/demultiplex/splitting/df2.txt")
        log:
            "logs/sample_tag_id2.log"
        shell:
            "paste -d '\t' {input.cb} {input.id} > {output.df2}"

    ## sorting both dataframes
    rule sample_tag_sort:
        input:
            df1="data/demultiplex/splitting/df1.txt",
            df2="data/demultiplex/splitting/df2.txt"
        output:
            df1_sorted=temp("data/demultiplex/splitting/df1_sorted.txt"),
            df2_sorted=temp("data/demultiplex/splitting/df2_sorted.txt")
        log:
            "logs/sample_tag_sort.log"
        threads: 8
        run:
            shell(
                """
                sort -t $'\t' -k1,1 {input.df2} -o {output.df2_sorted}
                sort -t $'\t' -k1,1 {input.df1} -o {output.df1_sorted}
                """
            )

    ## joining both dataframes --> dataframe: ceullar barcode + sample tag + ID
    rule sample_tag_join:
        input:
            df1_sorted="data/demultiplex/splitting/df1_sorted.txt",
            df2_sorted="data/demultiplex/splitting/df2_sorted.txt"
        output:
            df=temp("data/demultiplex/splitting/df.txt")
        log:
            "logs/sample_tag_join.log"
        shell:
            "join -t $'\t' -1 1 -2 1 {input.df2_sorted} {input.df1_sorted} > {output.df}"

    ## creating as many dataframes as sample tags; each dataframe contains IDs that belong to one specific sample tag
    rule sample_tag_df:
        input:
            "data/demultiplex/splitting/df.txt"
        output:
            directory("data/demultiplex/splitting/df")
        run:
            shell(
                """
                mkdir data/demultiplex/splitting/df
                awk -F' ' '{{ print $2 >> ("data/demultiplex/splitting/df/" $3 "_df.txt") }}' {input}
                """
            )

    ## splitting the original fastq file into as many fastqs as sample tags according to the ID dataframes created in the previous rule
    rule fastq_split:
        input: 
            fastq=config['raw_data_fastq_list'][1],
            df="data/demultiplex/splitting/df"
        output:
            directory("data/demultiplex/splitting/fastq")
        params: 
            threads_number = config['threads_number']
        run:
            shell(
                """
                mkdir data/demultiplex/splitting/fastq
                ls {input.df}/*_df.txt | xargs -I {{}} -P {params.threads_number} sh -c 'file="{{}}"; basefile=$(basename "$file"); seqtk subseq {input.fastq} "$file" | gzip > data/demultiplex/splitting/fastq/"${{basefile%_*}}.R2.fastq.gz"'
                """
            )

## ------------------------------------------------------------------------------------ ##
## amplicon-based data
## ------------------------------------------------------------------------------------ ##
if not config['wta']:

    ## ------------------------------------------------------------------------------------ ##
    ## allele-typing
    ## ------------------------------------------------------------------------------------ ##
    ## changing cDNA headers
    rule change_headers: 
        input:
            original_fasta=config['amplicon_cDNA_fasta']
        output:
            updated_cDNA="data/allele_typing/cDNA.fasta"
        shell:
            "cat {input} | sed 's/location.*AMPLICON//' > {output}"

    ## trimming sequences on which allele-typing will be performed to 15bp
    rule extract_sequences: 
        input:
            cdna_fasta="data/allele_typing/cDNA.fasta"
        output:
            short_seq="data/allele_typing/short_seq.fasta"
        params:
            genes=config['genes_to_be_allele_typed']
        run:
            shell(
                """
                fasta_file={input.cdna_fasta}
                gene_names=({params.genes})
                pattern=$(IFS='|'; echo "${{gene_names[*]}}")
                seq_data=$(awk -v pattern="$pattern" -v RS=">" -v FS="\\n" ' $1 ~ pattern {{ printf ">%s", $0 }} ' "$fasta_file")
                echo "$seq_data" | awk '/^>/ {{ header=$0; getline; sequence=substr($0, 1, 15); print header"\\n"sequence }}' > {output.short_seq}
                """
            )
    ## performing allele-typing
    rule allele_typing:
        input:
            fastq=lambda wildcards: "data/demultiplex/splitting/fastq" if config['multiplex'] else "data/raw",
            short_seq="data/allele_typing/short_seq.fasta",
            script="scripts/3_allele_typing.sh"
        output:
            directory("data/allele_typing/alleles")
        params:
            threads_number=config['threads_number']
        run:
            shell(
                """
                mkdir data/allele_typing/alleles
                chmod +x {input.script}
                ./{input.script} {input.fastq} {input.short_seq} {params.threads_number}
                """
            )

    ## ------------------------------------------------------------------------------------ ##
    ## quantification
    ## ------------------------------------------------------------------------------------ ##
    ## adding the typed-allele sequences into our final cDNA file
    rule final_cDNA_fasta:
        input:
            alleles_dir="data/allele_typing/alleles",
            genes_cDNA="data/allele_typing/cDNA.fasta",
            script="scripts/4_final_cDNA_amplicon.sh"
        output:
            cDNA="data/quant/cDNA.fasta",
            lookup_table="data/quant/lookup_table_HLA.csv"
        params:
            genes=config['genes_to_be_allele_typed']
        run:
            shell(
                """
                chmod +x {input.script}
                ./{input.script} {input.alleles_dir} {params.genes} {input.genes_cDNA} {output.cDNA} {output.lookup_table}
                """
            )

    ## creating a transcript to gene txt file
    rule create_t2g: 
        input:
            cDNA="data/quant/cDNA.fasta"
        output:
            t2g="data/quant/t2g.txt"
        run:
            shell(
                """
                file1=$(grep "^>" {input.cDNA} | sed 's/^.//')
                file2=$(cat {input.cDNA} | sed 's/|.*//' | grep "^>" | sed 's/^.//')
                paste <(echo "$file1") <(echo "$file2") > {output.t2g}
                """
            )
    
    ## creating cDNA index
    rule create_index:
        input:
            cDNA="data/quant/cDNA.fasta"
        output:
            index="data/quant/index.idx"
        shell:
            "kallisto index -i {output.index} {input.cDNA}"
    
    ## quantification using kallisto
    rule quantification:
        input:
            index="data/quant/index.idx",
            t2g="data/quant/t2g.txt",
            raw_data_fastq_list=config['raw_data_fastq_list']
        output:
            directory("data/quant/kb_output_amplicon")
        params:
            tech=config['sequencing_technology']
        log:
            "logs/quantification.log"
        resources:
            mem_gb=8,
            cores=8
        shell:
            """  
            kb count -i {input.index} -g {input.t2g} -x {params.tech} -o {output} --mm --verbose {input.raw_data_fastq_list} > {log} 2>&1
            mv data/quant/lookup_table_HLA.csv data/quant/kb_output_amplicon/counts_unfiltered/
            """

## ------------------------------------------------------------------------------------ ##
## whole transcriptome data
## ------------------------------------------------------------------------------------ ##
if config['wta']:

    ## ------------------------------------------------------------------------------------ ##
    ## allele-typing
    ## ------------------------------------------------------------------------------------ ##
    ## unzipping reference and gtf files
    rule gunzip:
        input:
            fasta=config['reference_genome_fasta'],
            gtf=config['reference_genome_gtf']
        output:
            fasta="data/meta/dna.primary_assembly.fa",
            gtf="data/meta/gtf.gtf"
        run:
            shell(
                """
                gunzip -c {input.fasta} > {output.fasta}
                gunzip -c {input.gtf} > {output.gtf}
                """
            )

    ## STAR genome generation
    rule genome_generation: 
        input:
            fasta="data/meta/dna.primary_assembly.fa",
            gtf="data/meta/gtf.gtf"
        output:
            genome_dir=directory("data/allele_typing/GenomeDir")
        params: 
            threads_number=config['threads_number']
        shell:
            "STAR --runMode genomeGenerate --genomeDir {output.genome_dir} --runThreadN {params.threads_number} --genomeSAindexNbases 14 --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf}"

    ## STAR alignment
    rule read_alignment:
        input:
            fastq=lambda wildcards: "data/demultiplex/splitting/fastq/*"[:2] if config['multiplex'] else config['raw_data_fastq_list'][:2],
            genome_dir="data/allele_typing/GenomeDir"
        output:
            star=temp("data/allele_typing/STAR_alignment/Aligned.out.bam")
        params:
            threads_number=config['threads_number']
        shell:
            "STAR --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outSAMtype BAM Unsorted --runThreadN {params.threads_number} --genomeDir {input.genome_dir} --readFilesCommand zcat --readFilesIn {input.fastq} --outFileNamePrefix data/allele_typing/STAR_alignment/"

    ## sorting BAM file
    rule bam_sort: 
        input:
            bam="data/allele_typing/STAR_alignment/Aligned.out.bam"
        output:
            sorted_bam="data/allele_typing/STAR_alignment/Aligned.out.sorted.bam"
        resources:
            mem_gb=8,
            cores=8
        shell:
            "samtools sort {input.bam} -o {output.sorted_bam}"

    ## indexing BAM file
    rule bam_index:
        input:
            sorted_bam="data/allele_typing/STAR_alignment/Aligned.out.sorted.bam"
        output:
            bam_index="data/allele_typing/STAR_alignment/Aligned.out.sorted.bam.bai"
        resources:
            mem_gb=4,
            cores=4
        shell:
            "samtools index {input.sorted_bam}"

    ## extracting reads mapping to chromosome 6
    rule reads_extract:
        input:
            bam="data/allele_typing/STAR_alignment/Aligned.out.sorted.bam",
            bam_index="data/allele_typing/STAR_alignment/Aligned.out.sorted.bam.bai"
        output:
            arcasHLA=directory("data/allele_typing/output")
        params:
            is_single=config['single_end_samples'],
            threads_number=config['threads_number']
        run:
            shell("arcasHLA reference --version 3.9.0 -v")
            if str(params.is_single).lower() == 'true':
                shell("arcasHLA extract {input.bam} -o {output.arcasHLA} -t {params.threads_number} -v --single")
            if str(params.is_single).lower() == 'false':
                shell("arcasHLA extract {input.bam} -o {output.arcasHLA} -t {params.threads_number} -v")

    ## performing allele-typing
    rule allele_typing:
        input:
            chr6_fastq="data/allele_typing/output/"
        output:
            arcasHLA="data/allele_typing/alleles/Aligned.genotype.json"
        params:
            genes_to_be_allele_typed=config['genes_to_be_allele_typed'],
            is_single=config['single_end_samples'],
            threads_number=config['threads_number']
        run:
            if str(params.is_single).lower() == 'true':
                shell(
                    """
                    genes=$(echo {params.genes_to_be_allele_typed} | sed 's/HLA-//g; s/ /,/g')
                    arcasHLA genotype {input.chr6_fastq}/*.fq.gz -o data/allele_typing/alleles/ -g "$genes" -t {params.threads_number} -v --single
                    """
                )
            if str(params.is_single).lower() == 'false':
                shell(
                    """
                    genes=$(echo {params.genes_to_be_allele_typed} | sed 's/HLA-//g; s/ /,/g')
                    arcasHLA genotype {input.chr6_fastq}/*.fq.gz -o data/allele_typing/alleles/ -g "$genes" -t {params.threads_number} -v
                    """
                )

    ## changing allele headers
    rule allele_headers:
        input:
            genotype="data/allele_typing/alleles/Aligned.genotype.json",
            allelelist="data/meta/Allelelist.txt",
            cdna="data/meta/hla_nuc.fasta",
            script="scripts/5_allele_headers.sh"
        output:
            allele_headers=temp("data/matching_headers.txt")
        run:
            shell(
                """
                chmod +x {input.script}
                ./{input.script} {input.genotype} {input.allelelist} {input.cdna} {output.allele_headers}
                """
            )

    ## creating a cDNA file from the typed-allele sequences
    rule allele_fasta:
        input:
            allele_headers="data/matching_headers.txt",
            cdna="data/meta/hla_nuc.fasta"
        output:
            allele_cdna="data/allele_typing/allele_cDNA.fasta",
            lookup_table="data/allele_typing/lookup_table_HLA.csv"
        shell:
            """
            awk 'NR==FNR{{a[">"$0];next}} /^>/{{f=$0 in a}} f' {input.allele_headers} {input.cdna} > {output.allele_cdna}
            echo "Allele,Gene,Function" > {output.lookup_table} && cat {output.allele_cdna} | awk '/^>/ {{ header=substr($0,2); split(header, fields, " "); gene = substr(fields[2], 1, index(fields[2], "*") - 1); print fields[2] "," "HLA-" gene}}' | awk -F, 'BEGIN {{OFS=","}} NR>0 {{if ($2=="HLA-A" || $2=="HLA-B" || $2=="HLA-C") $3="HLA_class_I"; else if ($2=="HLA-DRA" || $2=="HLA-DRB1" || $2=="HLA-DRB3" || $2=="HLA-DRB4" || $2=="HLA-DRB5" || $2=="HLA-DQA1" || $2=="HLA-DQA2" || $2=="HLA-DQB1" || $2=="HLA-DQB2" || $2=="HLA-DPA1" || $2=="HLA-DPA2" || $2=="HLA-DPB1" || $2=="HLA-DPB2" || $2=="HLA-DMA" || $2=="HLA-DMB") $3="HLA_class_II"; else if ($2=="HLA-E" || $2=="HLA-F" || $2=="HLA-G") $3="HLA_non_classical"}} 1' >> {output.lookup_table}
            """

    ## ------------------------------------------------------------------------------------ ##
    ## quantification
    ## ------------------------------------------------------------------------------------ ##
    ## creating a transcript to gene txt file & cDNA index & final cDNA
    rule ref_cDNA_fasta:
        input:
            fasta=config['reference_genome_fasta'],
            gtf=config['reference_genome_gtf']
        output:
            index="data/quant/meta/ref_index.idx",
            t2g="data/quant/meta/ref_t2g.txt",
            fasta=temp("data/quant/genes_cDNA.fa")
        resources:
            mem_gb=8,
            cores=8
        shell:
            "kb ref -i {output.index} -g {output.t2g} -f1 {output.fasta} {input.fasta} {input.gtf}"

    rule final_cDNA_fasta:
        input:
            allele_cDNA="data/allele_typing/allele_cDNA.fasta",
            genes_cDNA="data/quant/genes_cDNA.fa",
            script="scripts/6_final_cDNA_wta.sh"
        output:
            cDNA="data/quant/cDNA.fasta",
        params:
            genes=config['genes_to_be_allele_typed']
        resources:
            cores=4
        run:
            shell(
                """
                rm -r data/quant/meta
                genes=$(echo {params.genes} | sed 's/ /|/g')
                chmod +x {input.script}
                ./{input.script} {input.allele_cDNA} {input.genes_cDNA} "$genes" {output.cDNA}
                """
            )

    rule create_t2g: 
        input:
            cDNA="data/quant/cDNA.fasta"
        output:
            t2g="data/quant/t2g.txt"
        shell:
            """
            grep "^>" {input.cDNA} | sed 's/^.//' | tr " " "\t" > {output.t2g}
            """

    rule create_index:
        input:
            cDNA="data/quant/cDNA.fasta"
        output:
            index="data/quant/index.idx"
        resources:
            cores=4
        shell:
            "kallisto index -i {output.index} {input.cDNA}"

    ## quantification using kallisto
    rule quantification:
        input:
            index="data/quant/index.idx",
            t2g="data/quant/t2g.txt",
            raw_data_fastq_list=config['raw_data_fastq_list']
        output:
            directory("data/quant/kb_output_WTA")
        params:
            tech=config['sequencing_technology']
        log:
            "logs/quantification.log"
        resources:
            mem_gb=8,
            cores=8
        shell:
            """  
            kb count -i {input.index} -g {input.t2g} -x {params.tech} -o {output} --mm --verbose {input.raw_data_fastq_list} > {log} 2>&1
            mv data/allele_typing/lookup_table_HLA.csv data/quant/kb_output_WTA/counts_unfiltered/
            """

## ------------------------------------------------------------------------------------ ##
##                            Success and failure messages
## ------------------------------------------------------------------------------------ ##
onsuccess:
	print("Success! The Snakemake workflow is completed.")

onerror:
	print("Error! The Snakemake workflow aborted.")
