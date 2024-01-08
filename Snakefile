configfile: 'config.yaml'

rule all:
    input:
        "data/quant/kb_output_amplicon" if not config['wta'] else "data/quant/kb_output_WTA"

if config['multiplex']:
    
    # demultiplexing
    rule sample_tag_prep:
        input: 
            fasta=config['sample_tag_seqs']
        output:
            tsv="data/demultiplex/meta/sample_tag.tsv"
        log:
            "logs/sample_tag_prep.log"
        script:
            "scripts/1_sample_tag_prep.R"

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

    rule sample_tag_ident:
        input:
            dir="data/demultiplex/quant"
        output:
            df1=temp("data/demultiplex/splitting/df1.txt")
        log:
            "logs/sample_tag_ident.log"
        script:
            "scripts/2_sample_tag_ident.R"

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

if not config['wta']:

    # amplicon allele-typing
    rule change_headers: 
        input:
            original_fasta=config['amplicon_cDNA_fasta']
        output:
            updated_cDNA="data/allele_typing/cDNA.fasta"
        shell:
            "cat {input} | sed 's/location.*AMPLICON//' > {output}"

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

    # amplicon quantification
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

    rule create_index:
        input:
            cDNA="data/quant/cDNA.fasta"
        output:
            index="data/quant/index.idx"
        shell:
            "kallisto index -i {output.index} {input.cDNA}"

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
            "kb count -i {input.index} -g {input.t2g} -x {params.tech} -o {output} --mm --verbose {input.raw_data_fastq_list} > {log} 2>&1"

if config['wta']:

    # WTA allele-typing
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

    rule genome_generation: 
        input:
            fasta="data/meta/dna.primary_assembly.fa",
            gtf="data/meta/gtf.gtf"
        output:
            genome_dir=temp("data/allele_typing/GenomeDir")
        params: 
            threads_number=config['threads_number']
        shell:
            "STAR --runMode genomeGenerate --genomeDir {output.genome_dir} --runThreadN {params.threads_number} --genomeSAindexNbases 14 --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf}"

    rule read_alignment:
        input:
            fastq=lambda wildcards: "data/demultiplex/splitting/fastq/*" if config['multiplex'] else config['raw_data_fastq_list'],
            genome_dir="data/allele_typing/GenomeDir"
        output:
            star="data/allele_typing/STAR_alignment/"
        params:
            threads_number=config['threads_number']
        shell:
            "STAR --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --runThreadN {params.threads_number} --genomeDir {input.genome_dir} --readFilesCommand zcat --readFilesIn {input.fastq} --outFileNamePrefix {output.star}"

    rule sam_to_bam:
        input:
            sam="data/allele_typing/STAR_alignment/Aligned.out.sam"
        output:
            bam=temp("data/allele_typing/STAR_alignment/Aligned.out.bam")
        resources:
            mem_gb=8,
            cores=8
        run:
            shell(
                """
                samtools view -b -S -o {output.bam} {input.sam}
                rm data/allele_typing/STAR_alignment/Aligned.out.sam
                """
            )

    rule bam_sort: 
        input:
            bam="data/allele_typing/STAR_alignment/Aligned.out.bam"
        output:
            sorted_bam=temp("data/allele_typing/STAR_alignment/Aligned.out.sorted.bam")
        resources:
            mem_gb=8,
            cores=8
        shell:
            "samtools sort {input.bam} -o {output.sorted_bam}"

    rule bam_index:
        input:
            sorted_bam="data/allele_typing/STAR_alignment/Aligned.out.sorted.bam"
        output:
            bam_index=temp("data/allele_typing/STAR_alignment/Aligned.out.sorted.bam.bai")
        resources:
            mem_gb=4,
            cores=4
        shell:
            "samtools index {input.sorted_bam}"

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
            shell(
                """
                arcasHLA reference --version 3.9.0 -v
                arcasHLA extract {input.bam} -o {output.arcasHLA} -t {params.threads_number} { "--single" if params.is_single else "" } -v
                """
            )

    rule allele_typing:
        input:
            chr6_fastq="data/allele_typing/output/"
        output:
            arcasHLA=directory("data/allele_typing/alleles")
        params:
            genes_to_be_allele_typed=config['genes_to_be_allele_typed'],
            is_single=config['single_end_samples'],
            threads_number=config['threads_number']
        run:
            shell(
                """
                genes=$(echo {params.genes_to_be_allele_typed} | sed s/HLA-//g; s/ /,/g')
                arcasHLA genotype {input.chr6_fastq}/*.fq.gz -o {output.arcasHLA} -g "$genes" -t {params.threads_number} { "--single" if params.is_single else "" } -v
                """
            )

    rule allele_headers:
        input:
            genotype="data/allele_typing/alleles",
            allelelist="data/meta/Allelelist.txt",
            cdna="data/meta/hla_nuc.fasta",
            script="scripts/5_allele_headers.sh"
        output:
            allele_headers=temp("data/matching_headers.txt")
        run:
            shell(
                """
                chmod +x {input.script}
                ./{input.script} {input.genotype}/Aligned.genotype.json {input.allelelist} {input.cdna} {output.allele_headers}
                """
            )

    rule allele_fasta:
        input:
            allele_headers="data/matching_headers.txt",
            cdna="data/meta/hla_nuc.fasta"
        output:
            allele_cdna="data/allele_typing/allele_cDNA.fasta"
        shell:
            """
            awk 'NR==FNR{{a[">"$0];next}} /^>/{{f=$0 in a}} f' {input.allele_headers} {input.cdna} > {output.allele_cdna}    
            """

    # WTA quantification
    rule ref_cDNA_fasta:
        input:
            fasta=config['reference_genome_fasta'],
            gtf=config['reference_genome_gtf']
        output:
            index="data/quant/meta/ref_index.idx",
            t2g="data/quant/meta/ref_t2g.txt",
            fasta="data/quant/genes_cDNA.fa"
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
            "kb count -i {input.index} -g {input.t2g} -x {params.tech} -o {output} --mm --verbose {input.raw_data_fastq_list} > {log} 2>&1"

## ------------------------------------------------------------------------------------ ##
##                            Success and failure messages
## ------------------------------------------------------------------------------------ ##
onsuccess:
	print("Success! The Snakemake workflow is completed.")

onerror:
	print("Error! The Snakemake workflow aborted.")