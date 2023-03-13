cd  /Users/admin/Analysis/06_201903_Triplex_Integrated/NR4A1_peak/Cancer-ATAC-ATLAS/peak_matrix/

# get only peaks with ID 

awk '$1=="PRAD_72732" || $1=="sample" {print}' TCGA_ATAC_peak_Log2Counts_dedup_sample  >PRAD_72732_signal_dedup.txt  

awk '$1=="PRAD_72733" || $1=="sample" {print}' TCGA_ATAC_peak_Log2Counts_dedup_sample  >PRAD_72733_signal_dedup.txt  




######################################################################################
