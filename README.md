# Metagenomics_Snakemake_Pipeline
Metagenomics snakemake pipeline to find Relative Species Abundance. 
 * The goal of this pipeline is to join all the following steps into one cohesive python script.
   
Development of pipeline in progress.  
  
Completed steps:  
 1)‎ Deinterleave with BBmap  
 2) Quality Check and Trimming with fastp   
 3) Assembly with MetaSPAdes
 4) Binning with MaxBin   
 5) Alignment with bwa index   
 6) Alignment with bwa alignment & conversion to bam format. 

Next steps:  
 7) Binning with MetaBat   
 8) Binning with DAS_Tool   
 9) Taxonomic Identification with GTDB-tk   
 10) Relative Species Abundance with CheckM   

  
