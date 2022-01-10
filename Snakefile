

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
        "plots/quals.svg"


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
