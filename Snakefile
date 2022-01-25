
#This Pipeline is meant to find relative species abundance.

#Current version with does not include filter out algal reads step.
#this may be important for genomics classification.
#pipeline may have to be run twice, once with algal reads and once without
#in order to find relative species abundance as well as individual species types.

#Helpful snakemake arguments if unstable cluster: --restart-times 10 --rerun-incomplete --keep-going

# Currently sample names are hard-coded.
# All samples should be in '.fastq.gz' format. The sample names without this suffix 
# should be present in the directory titled 'data'. Nothing else should be in this directory.

SAMPLES = ["ChrtobMetaP3_0mM_FD", "ChrtobMetaP5_0mM_FD"]

#Set threads based on availability.
THREADS = 14

#If using slurm / sbatch: gtdbtk requires extensive amounts of memory. Between 300-500GB is reccomended per sample.
#Do not run this locally! OOM-kill at best.

#Here are all of the expected outputs. 
rule all:
    input:
        expand("OD_summary_deinterleave/{sample}_deinterleave_summary.txt", sample = SAMPLES),
        expand("OD_summary_fastp/{sample}_fastp_summary.txt", sample = SAMPLES),
        expand("OD_assembled/{sample}_assembled/corrected/{sample}.out.R2.fq.00.0_0.cor.fastq.gz", sample = SAMPLES),
        expand("OD_maxbin/{sample}_log" , sample = SAMPLES),
        expand("OD_bwa_index/{sample}.sa", sample = SAMPLES),
        expand("OD_aligned/{sample}_aligned.bam", sample = SAMPLES),
        expand("OD_metabat/{sample}_log", sample = SAMPLES),
        expand("OD_scaffolds2bin/{sample}_metabat.scaffolds2bin.tsv", sample = SAMPLES),
        expand("OD_dastool/{sample}_DASTool_summary.txt", sample = SAMPLES),
        "OD_gtdbtk/gtdbtk_database_download_log",
        expand("OD_gtdbtk/{sample}_log", sample = SAMPLES),
        expand("OD_aligned/{sample}_aligned.bam.bai", sample = SAMPLES),
        expand("OD_checkm_coverage/{sample}_log", sample = SAMPLES),
        expand("OD_checkm_profile/{sample}_profile.txt", sample = SAMPLES)


# Deinterleave fastq files.
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
        module load bbmap
        reformat.sh \
        in={input} \
        out1={output[0]} \
        out2={output[1]} \
        2>> {output[2]}
        """

#Trim, cut adapters, quality filter.
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

#Assembly of Genomes
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

#Bin genomes
rule maxbin:
    "Run Maxbin"
    input:
        "OD_assembled/{sample}_assembled/scaffolds.fasta",
        "OD_fastp/{sample}.out.R1.fq.gz",
        "OD_fastp/{sample}.out.R2.fq.gz"
    output:
        "OD_maxbin/{sample}_log",
        "OD_maxbin/{sample}_binned.log"
    threads: THREADS
    shell:
        """
        run_MaxBin.pl \
        -contig {input[0]} \
        -reads {input[1]} \
        -reads2 {input[2]} \
        -out OD_maxbin/{wildcards.sample}_binned \
        -thread {threads} 2> {output[0]}
        """

#Make index for BWA alignment.
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
#BWA Alignment step. This is required for metabat.
rule bwa_alignment:
    "Run bwa alignment"
    input:
        "OD_bwa_index/{sample}.sa",
        "OD_assembled/{sample}_assembled/corrected/{sample}.out.R1.fq.00.0_0.cor.fastq.gz",
        "OD_assembled/{sample}_assembled/corrected/{sample}.out.R2.fq.00.0_0.cor.fastq.gz"
    output:
        "OD_aligned/{sample}_aligned.bam"
    threads: THREADS
    shell:
        """
        bwa mem -t {threads} \
        OD_bwa_index/{wildcards.sample} \
        {input[1]} \
        {input[2]} \
        | samtools view -S -b | samtools sort > {output} 
        """
#This is another binning tool.
rule Metabat:
    "Run Metabat"
    input:
        "OD_assembled/{sample}_assembled/scaffolds.fasta",
        "OD_aligned/{sample}_aligned.bam"
    output:
        "OD_metabat/{sample}_log",
        "{sample}_depth.txt"
    threads: 12
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output[1]} {input[1]}
        metabat2 -t {threads} -i {input[0]} -a {output[1]} -o OD_metabat/{wildcards.sample} 2> {output[0]}
        """
#Restructure database to more easily perform next step.
rule restructure_bins:
    "Restructure bins"
    input:
         "OD_maxbin/{sample}_log",
         "OD_metabat/{sample}_log",
         "{sample}_depth.txt"
    output:
        directory("OD_maxbin/{sample}_bins/"),
        directory("OD_metabat/{sample}_bins/")
    threads: 1
    shell:
        """
        mkdir OD_metabat/{wildcards.sample}_bins
        mv OD_metabat/{wildcards.sample}.* OD_metabat/{wildcards.sample}_bins/
        mkdir OD_maxbin/{wildcards.sample}_bins
        mv OD_maxbin/{wildcards.sample}_binned.0* OD_maxbin/{wildcards.sample}_bins/
        mv {wildcards.sample}_depth.txt ./OD_metabat/
        """

#Prepare for dastool.
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

#Dastool takes the previously found bins from maxbin and metabat, and chooses the best / dereplicates them.
rule dastool:
    "Run DAS_Tool"
    input:
        "OD_scaffolds2bin/{sample}_maxbin.scaffolds2bin.tsv",
        "OD_scaffolds2bin/{sample}_metabat.scaffolds2bin.tsv",
        "OD_assembled/{sample}_assembled/scaffolds.fasta"
    output:
        "OD_dastool/{sample}_log",
        "OD_dastool/{sample}_DASTool_summary.txt",
        directory("OD_dastool/{sample}_DASTool_bins/")
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

#database installation for gtdbtk
#this rule will fail if database is already installed. 
#Comment this out if database is already installed. Also comment out its input in Rule ALL.
#database is large. 50 or more GB.
rule gtdbtk_database:
    "Download gtdbtk database"
    output:
        "OD_gtdbtk/gtdbtk_database_download_log"
    threads: THREADS
    shell:
        """
        download-db.sh 2> {output}
        """

#highly memory intensive, most common reason for oom-kill.
#This classifies genomes based on the bins. Should determine what each species is.
rule gtdbtk:
    "Run gtdbtk"
    input:
        "OD_dastool/{sample}_DASTool_bins/",
        "OD_gtdbtk/gtdbtk_database_download_log"
    output:
        "OD_gtdbtk/{sample}_log"
    threads: THREADS
    shell:
        """
        gtdbtk classify_wf \
        --genome_dir {input[0]} \
        --out_dir OD_gtdbtk/{wildcards.sample}/ \
        --cpus {threads} \
        --extension fa 2> {output}
        """

#required for checkm
rule index_bam:
    "Run index bam file before checkM"
    input:
        "OD_aligned/{sample}_aligned.bam"
    output:
        "OD_aligned/{sample}_aligned.bam.bai",
        "OD_aligned/{sample}_index_log"
    threads: 4
    shell:
        """
        samtools index {input} {output[0]} 2> {output[1]}
        """


#mamba install -c bioconda/label/cf201901 checkm-genome
#This checkm install command will also install reference data / set root for checkm.
#other checkm installations may not include reference data.

#Checkm coverage is needed to run checkm profile
rule checkm_coverage:
    "Run checkm coverage"
    input:
        "OD_dastool/{sample}_DASTool_bins/",
        "OD_aligned/{sample}_aligned.bam",
        "OD_aligned/{sample}_aligned.bam.bai"
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

#Checkm profile determines the relative species abundance
#Current version of this package includes a bug:
#To fix bug, in checkm 'profile.py' file:
#change "import checkm.prettytable" to "import checkm.prettytable as prettytable"
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
