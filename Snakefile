
SAMPLES = ["ChrtobMetaP3_0mM_FD", "ChrtobMetaP5_0mM_FD", "ChrtobMetaP5_4mM_FD"]
THREADS = 12

# rule all:
#     input:

#         OD_summary_fastp # this checks to make sure rule 1 output was successful before moving to rule 2
        # OD_summary_metaspades, # and so on...
        # OD_summary_bwa_alignment,
        # OD_summary_gtdbtk,
        # OD_summary_checkM
       # expand("{sample}.R{read_no}.fq.gz.out", sample=SAMPLES, read_no=['1', '2'])


# OD_deinterleave          = OUT + "/01_deinterleave"
# OD_summary_deinterleave  = OUT + "/02_summary_deinterleave"
# OD_fastp                 = OUT + "/03_fastp"
# OD_summary_fastp         = OUT + "/04_summary_fastp"

rule all:
    input:
        expand("OD_summary_deinterleave/{sample}_deinterleave_summary.txt", sample = SAMPLES),
        expand("OD_summary_fastp/{sample}_fastp_summary.txt", sample = SAMPLES),
        expand("OD_assembled/{sample}_assembled/corrected/{sample}.out.R2.fq.00.0_0.cor.fastq.gz", sample = SAMPLES),
        expand("OD_maxbin/{sample}_log" ,sample = SAMPLES),
        expand("OD_bwa_index/{sample}.sa", sample = SAMPLES),
        expand("OD_aligned/{sample}_aligned.bam", sample = SAMPLES),
        expand("OD_metabat/{sample}_log", sample = SAMPLES),
        expand("OD_scaffolds2bin/{sample}_metabat.scaffolds2bin.tsv", sample = SAMPLES),
        expand("OD_dastool/{sample}_DASTool_summary.txt", sample = SAMPLES),
        "OD_gtdbtk/gtdbtk_database_download_log",
        expand("OD_gtdbtk/{sample}_log", sample = SAMPLES),
        expand("OD_aligned/{sample}_index_log", sample = SAMPLES),
        expand("OD_checkm_coverage/{sample}_log", sample = SAMPLES),
        expand("OD_checkm_profile/{sample}_profile.txt", sample = SAMPLES)



#load modules:
#"module load bbmap"

rule deinterleave:
    "Run BBmap reformat.sh"
    input:
        "data/{sample}.fastq.gz"
    output:
        "OD_deinterleave/{sample}_read1.fq",
        "OD_deinterleave/{sample}_read2.fq",
        "OD_summary_deinterleave/{sample}_deinterleave_summary.txt"
    shell:
        """
        reformat.sh \
        in={input} \
        out1={output[0]} \
        out2={output[1]} \
        2>> {output[2]}
        """

rule fastp:
    "Run fastp"
    input:
        "OD_deinterleave/{sample}_read1.fq",
        "OD_deinterleave/{sample}_read2.fq"
    output:
        "OD_fastp/{sample}.report.json",
        "OD_fastp/{sample}.report.html",
        "OD_fastp/{sample}.out.R1.fq.gz",
        "OD_fastp/{sample}.out.R2.fq.gz",
        "OD_summary_fastp/{sample}_fastp_summary.txt"
    threads: THREADS
    shell:
        """
        fastp \
        -i {input[0]} \
        -I {input[1]} \
        -c \
        -w {threads} \
        --cut_front \
        --cut_front_window_size=1 \
        --cut_front_mean_quality=3 \
        --cut_tail \
        --cut_tail_window_size=1 \
        --cut_tail_mean_quality=3 \
        --cut_right \
        --cut_right_window_size=15 \
        --cut_right_mean_quality=15 \
        -l 35 \
        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -j {output[0]} \
        -h {output[1]} \
        -o {output[2]} \
        -O {output[3]} \
        1>> {output[4]}
        """

rule MetaSPAdes:
    "Run Metaspades Assembly"
    input:
        "OD_fastp/{sample}.out.R1.fq.gz",
        "OD_fastp/{sample}.out.R2.fq.gz"
    output:
        "OD_assembled/{sample}_assembled/scaffolds.fasta",
        "OD_assembled/{sample}_assembled/corrected/{sample}.out.R1.fq.00.0_0.cor.fastq.gz",
        "OD_assembled/{sample}_assembled/corrected/{sample}.out.R2.fq.00.0_0.cor.fastq.gz"
    threads: THREADS
    shell:
        """
       metaspades.py \
        -o "OD_assembled/{wildcards.sample}_assembled/" \
        -1 {input[0]} \
        -2 {input[1]} \
        -k 21,33,55,77,99,127 \
        --threads {threads}
        """

#note for maxbin. installing this led to package changes. 
rule maxbin:
    "Run Maxbin"
    input:
        "OD_assembled/{sample}_assembled/scaffolds.fasta",
        "OD_fastp/{sample}.out.R1.fq.gz",
        "OD_fastp/{sample}.out.R2.fq.gz"
    output:
        "OD_maxbin/{sample}_log",
        "OD_maxbin/{sample}_bins"
    threads: THREADS
    shell:
        """
        run_MaxBin.pl \
        -contig {input[0]} \
        -reads {input[1]} \
        -reads2 {input[2]} \
        -out OD_maxbin/{wildcards.sample}_binned/ \
        -thread {threads} 2> {output}
        sleep 10s
        mkdir OD_maxbin/{wildcards.sample}_bins
        sleep 10s
        mv OD_maxbin/{wildcards.sample}_binned.0* OD_maxbin/{wildcards.sample}_bins/
        sleep 5s
        """

rule bwa_index:
    "Run bwa indexer"
    input:
        "OD_assembled/{sample}_assembled/scaffolds.fasta"
    output:
        "OD_bwa_index/{sample}.sa"
    threads: 1
    shell:
        """
        bwa index -p OD_bwa_index/{wildcards.sample} {input} 
        """

rule bwa_alignment:
    "Run bwa alignment"
    input:
        "OD_assembled/{sample}_assembled/corrected/{sample}.out.R1.fq.00.0_0.cor.fastq.gz",
        "OD_assembled/{sample}_assembled/corrected/{sample}.out.R2.fq.00.0_0.cor.fastq.gz"
    output:
        "OD_aligned/{sample}_aligned.bam"
    threads: THREADS
    shell:
        """
        bwa mem -t {threads} \
        OD_bwa_index/{wildcards.sample} \
        {input[0]} \
        {input[1]} \
        | samtools view -S -b | samtools sort > {output} 
        """

rule Metabat:
    "Run Metabat"
    input:
        "OD_assembled/{sample}_assembled/scaffolds.fasta",
        "OD_aligned/{sample}_aligned.bam"
    output:
        "OD_metabat/{sample}_log",
        "OD_metabat/{sample}_bins"
    threads: THREADS
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {wildcards.sample}_depth.txt {input[1]}
        metabat2 -t={threads} -i {input[0]} -a {wildcards.sample}_depth.txt -o "OD_metabat/{wildcards.sample}" 2> {output}
        mkdir OD_metabat/{wildcards.sample}_bins
        mv OD_metabat/{wildcards.sample}.* OD_metabat/{wildcards.sample}_bins/
        """


rule fasta_scaffold2bin:
    "Run Fasta_to_Scaffolds2Bin.sh - this is required for dastool input"
    input:
        "OD_maxbin/{sample}_bins/",
        "OD_metabat/{sample}_bins/"
    output:
        "OD_scaffolds2bin/{sample}_maxbin.scaffolds2bin.tsv",
        "OD_scaffolds2bin/{sample}_metabat.scaffolds2bin.tsv"
    threads: THREADS
    shell:
        """
        Fasta_to_Scaffolds2Bin.sh -i {input[0]} -e fasta > {output[0]}
        Fasta_to_Scaffolds2Bin.sh -i {input[1]} -e fa > {output[1]}
        """


rule dastool:
    "Run DAS_Tool"
    input:
        "OD_scaffolds2bin/{sample}_maxbin.scaffolds2bin.tsv",
        "OD_scaffolds2bin/{sample}_metabat.scaffolds2bin.tsv",
        "OD_assembled/{sample}_assembled/scaffolds.fasta"
    output:
        "OD_dastool/{sample}_log",
        "OD_dastool/{sample}_DASTool_summary.txt",
        "OD_dastool/{sample}_DASTool_bins"
    threads: THREADS
    shell:
        """
        DAS_Tool -i {input[0]},{input[1]} \
        -l maxbin,metabat \
        -c {input[2]} \
        -o OD_dastool/{wildcards.sample} \
        -t {threads} \
        --search_engine diamond \
        --write_bins 1 2> {output[0]}
        """

rule gtdbtk_database:
    "Download gtdbtk database"
    output:
        "OD_gtdbtk/gtdbtk_database_download_log"
    threads: THREADS
    shell:
        """
        download-db.sh 2> {output}
        """


#--restart-times 30 --rerun-incomplete --keep-going


rule gtdbtk:
    "Run gtdbtk"
    input:
        "OD_dastool/{sample}_DASTool_bins/"
    output:
        "OD_gtdbtk/{sample}_log"
    threads: THREADS
    shell:
        """
        gtdbtk classify_wf \
        --genome_dir {input} \
        --out_dir OD_gtdbtk/{wildcards.sample}/ \
        --cpus {threads} \
        --extension fa 2> {output}
        """


rule index_bam:
    "Run index bam file before checkM"
    input:
        "OD_aligned/{sample}_aligned.bam"
    output:
        "OD_aligned/{sample}_index_log"
    threads: THREADS
    shell:
        """
        samtools index {input} {input}.bai 2> {output}
        """


# rule checkm_reference_data:
#     output:
#         "OD_checkm_reference_data/checkm_database_download_log"
#     threads: THREADS
#     shell:
#         """
#         wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz 2> {output}
#         mv checkm_data_2015_01_16.tar.gz OD_checkm_reference_data/
#         tar -xf OD_checkm_reference_data/checkm_data_2015_01_16.tar.gz
#         checkm data setRoot OD_checkm_reference_data
#         """

#mamba install -c bioconda/label/cf201901 checkm-genome

rule checkm_coverage:
    "Run checkm coverage"
    input:
        "OD_dastool/{sample}_DASTool_bins/",
        "OD_aligned/{sample}_aligned.bam"
    output:
        "OD_checkm_coverage/{sample}_log",
        "OD_checkm_coverage/{sample}_checkm_cov_out.tsv"
    threads: THREADS
    shell:
        """
        checkm coverage \
        -x fa \
        -t {threads} \
        {input[0]} \
        {output[1]} \
        {input[1]} \
        2> {output[0]}
        """


rule checkm_profile:
    "Run checkm profile"
    input:
        "OD_checkm_coverage/{sample}_checkm_cov_out.tsv"
    output:
        "OD_checkm_profile/{sample}_profile.txt"
    threads: THREADS
    shell:
        """
        checkm profile \
        -f {output} \
        {input}
        """






















