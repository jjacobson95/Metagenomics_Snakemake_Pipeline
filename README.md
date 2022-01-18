# Metagenomics_Snakemake_Pipeline
Metagenomics snakemake pipeline to find Relative Species Abundance. 
 * The goal of this pipeline is to join all the following steps into one cohesive python script.
   
Development of pipeline in progress.  
  
Completed steps:  
 1)‎ Deinterleave with **BBmap**  
 2) Quality Check and Trimming with **fastp**   
 3) Assembly with **MetaSPAdes** 
 4) Binning with **MaxBin**   
 5) Create **bwa index** files using 'bwa index'  
 6) Alignment with **bwa mem** & conversion to bam format using **samtools**  
 7) Binning with **MetaBat**   
 8) Preparation for **DAS_Tool** with **Fasta_to_Scaffolds2Bin.sh**. 
 9) Binning with **DAS_Tool**   
   
 Next Steps:  
 10) Taxonomic Identification with **GTDB-tk**   
   * This step is computationally intense and requires a substantial amount of memory.
 11) Relative Species Abundance with **CheckM**   

  
