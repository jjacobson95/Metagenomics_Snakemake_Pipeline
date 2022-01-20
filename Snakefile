
SAMPLES = ["ChrtobMetaP3_0mM_FD", "ChrtobMetaP5_0mM_FD", "ChrtobMetaP5_4mM_FD"]


# rule all:
#     input:
#         OD_summary_fastp # this checks to make sure rule 1 output was successful before moving to rule 2
        # OD_summary_metaspades, # and so on...
        # OD_summary_bwa_alignment,
        # OD_summary_gtdbtk,
        # OD_summary_checkM
       # expand("{sample}.R{read_no}.fq.gz.out", sample=SAMPLES, read_no=['1', '2'])

rule all:
    input:
        "data/"
    output:
         "OD_dastool/{wildcards.sample}_log"



#load modules:
#"module load bbmap"

rule deinterleave:
    "Run BBmap reformat.sh "
    input:
        "data/{sample}.fastq.gz"
    output:
        "OD_deinterleave/{sample}_read1.fq",
        "OD_deinterleave/{sample}_read2.fq"
    shell:
        """
        reformat.sh \
        in={input} \
        out1={output[0]} \
        out2={output[1]}
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
        "OD_fastp/{sample}.out.R2.fq.gz"
    shell:
        """
        fastp \
        -i {input[0]} \
        -I {input[1]} \
        -c \
        -w 8 \
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
        -O {output[3]} 
        """
 #this step is in progress. Check slurm file upon return.
rule MetaSPAdes:
    "Run Metaspades Assembly"
    input:
        "OD_fastp/{sample}.out.R1.fq.gz",
        "OD_fastp/{sample}.out.R2.fq.gz"
    output:
        "OD_assembled/{sample}_assembled/"
    shell:
        """
       metaspades.py \
        -o {output} \
        -1 {input[0]} \
        -2 {input[1]} \
        -k 21,33,55,77,99,127 \
        --threads 40
        """

#note for maxbin. installing this led to package changes. 
rule Maxbin:
    "Run Maxbin"
    input:
        "OD_assembled/{sample}_assembled/scaffolds.fasta",
        "OD_fastp/{sample}.out.R1.fq.gz",
        "OD_fastp/{sample}.out.R2.fq.gz"
    output:
        "OD_maxbin/{sample}_binned/",
        "OD_maxbin/{sample}_log"
    shell:
        """
        run_MaxBin.pl \
        -contig {input[0]} \
        -reads {input[1]} \
        -reads2 {input[2]} \
        -out {output[0]} \
        -thread 8 2> {output[1]}
        mkdir OD_maxbin/{wildcards.sample}_bins
        mv OD_maxbin/{wildcards.sample}_binned.0* OD_maxbin/{wildcards.sample}_bins/
        """

rule bwa_index:
    "Run bwa indexer"
    input:
        "OD_assembled/{sample}_assembled/scaffolds.fasta"
    output:
        "OD_bwa_index/{sample}"
    shell:
        """
        bwa index -p {output} {input} 
        """

rule bwa_alignment:
    "Run bwa alignment"
    input:
        "OD_assembled/{sample}_assembled/corrected/{sample}.out.R1.fq.00.0_0.cor.fastq.gz",
        "OD_assembled/{sample}_assembled/corrected/{sample}.out.R2.fq.00.0_0.cor.fastq.gz"
    output:
        "OD_aligned/{sample}_aligned.bam"
    shell:
        """
        bwa mem -t 20 \
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
        "OD_metabat/{sample}_log"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {wildcards.sample}_depth.txt {input[1]}
        metabat2 -i {input[0]} -a {wildcards.sample}_depth.txt -o "OD_metabat/{wildcards.sample}" 2> {output}
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
        "OD_dastool/{sample}_log"
    shell:
        """
        DAS_Tool -i {input[0]},{input[1]} \
        -l maxbin,metabat \
        -c {input[2]} \
        -o OD_dastool/{wildcards.sample} \
        -t 8 \
        --search_engine diamond \
        --write_bins 1 \
        2> {output}
        """

rule gtdbtk_database:
    "Download gtdbtk database"
    output:
        "OD_gtdbtk/gtdbtk_database_download_log"
    shell:
        """
        download-db.sh 2> {output}
        """


rule gtdbtk:
    "Run gtdbtk"
    input:
        "OD_dastool/{sample}_DASTool_bins/"
    output:
        "OD_gtdbtk/{sample}_log"
    shell:
        """
        gtdbtk classify_wf \
        --genome_dir {input} \
        --out_dir OD_gtdbtk/{wildcards.sample}/ \
        --cpus 40 \
        --extension fa 2> {output}
        """

















