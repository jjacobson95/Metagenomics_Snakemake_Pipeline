# Metagenomics_Snakemake_Pipeline
Metagenomics pipeline to find Relative Species Abundance. 
   
Pipeline Now Functional. View visualized workflow in metagenomics_dag.pdf
  
 Steps:  
 1)‎ Deinterleave with **BBmap**  
 2) Quality Check and Trimming with **fastp**   
 3) Assembly with **MetaSPAdes**  
 4) Binning with **MaxBin**   
 5) Create **bwa index** files using 'bwa index'  
 6) Alignment with **bwa mem** & conversion to bam format using **samtools**  
 7) Binning with **MetaBat**   
 8) Restructure File system 
 9) Preparation for **DAS_Tool** with **Fasta_to_Scaffolds2Bin.sh**  
 10) Binning with **DAS_Tool**   
 11) Download **GTDB-tk database** with **download-db.sh**  
 12) Taxonomic Identification with **GTDB-tk**   
 13) Index bam file with **samtools** for **CheckM**    
 14) Run **checkm coverage** to create coverage tsv in **CheckM**  
 15) Relative Species Abundance is found using **checkm profile** with **CheckM**   

  

See notes in Snakefile for more information.
