#!/bin/bash

#Run as: bash ./01_ChIRP-seq_analysis.sh

set -e

echo "----- ChIRP-seq PIPELINE (SE) -----"

project_folder="/01_Projects/lncRNA_project/ATAC-seq_HeLa_24h"
annotation="/01_Projects/genome/gencode.v27.annotation+chrR-Reversed_bin500.gtf"
genome="/01_Projects/genome/Bowtie2_hg38+chrR-rev_index/hg38+chrR"
fasta="/01_Projects/genome/Hg38_chrR-Reversed_bigOnly.fa"
blacklistedRegions="/01_Projects/genome/GRCh38.blacklist.bed"
chromSizes="/01_Projects/genome/hg38_chrR_U13369_bigChr.chrom.sizes"

echo  "$BLUE--------- Creatwe Bowtie2 "hg38+chrR" index ---------$RESET"

bowtie2-build $fasta $genome

echo "----- FastQC analysis -----"

cd $project_folder

mkdir -p raw_data results/MACS2 results/featureCount results/fastQC/ processed_data/aligned/ logs scripts 

fastqc $project_folder/raw_data/*.fastq.gz -o $project_folder/results/fastQC/ -t 9
multiqc $project_folder/results/fastQC/ -o $project_folder/results/fastQC/

echo  "$BLUE--------- Alignment to the refrence genome ---------$RESET"

cd $project_folder/raw_data

for file in $project_folder/raw_data/*.fastq.gz; do echo "Processing sample: " $file ; \
 bowtie2 --threads 16 --local -t --very-sensitive-local -x $genome -1 $file 2>> $project_folder/logs/Bowtie2_alignment_stats.txt | \
 samtools view -bS -f 2 -@ 8 - > $project_folder/processed_data/aligned/${file/%_1.fastq.gz/_PE_aligned.bam}; done

for file in $project_folder/aligned/*_aligned.bam; do echo "Processing sample: " $file; \
 samtools sort -@ 20 $file > ${file/.bam/.sorted.bam} ; done

for file in $project_folder/aligned/*.sorted.bam; \
 do echo "Processing sample: " $file; samtools index -b $file ;done

# Filter reads with low mapping quality 
for file in $project_folder/results/aligned/*.bam; do echo $file; alignmentSieve -b $file -o ${file/.bam/.filtered.bam}
 --minMappingQuality 20 --filterMetrics ${file/.bam/.log}; done


echo  "$BLUE--------- Peak calling ---------$RESET"

#Merged replicates - IP
for IP in $project_folder/results/aligned/IP_*.bam; do echo $IP; macs2 callpeak -t $IP \
	-c $project_folder/results/aligned/Input_*.bam -n H3K27ac_IP_merged -f BAM -g hs --nomodel --extsize 147 \
    --outdir $project_folder/results/MACS2 2> $project_folder/logs/H3K27ac_IP_merged_MACS2.log; done

#Merged replicates - Negative IP
for NegIP in $project_folder/results/aligned/NegIP_*.bam; do echo $NegIP; macs2 callpeak -t $NegIP \
	-c $project_folder/results/aligned/Input_*.bam -n H3K27ac_negIP_merged -f BAM -g hs --nomodel --extsize 147 \
    --outdir $project_folder/results/MACS2 2> $project_folder/logs/H3K27ac_negIP_merged_MACS2.log; done

echo  "$BLUE--------- Peak annotation ---------$RESET"

annotatePeaks.pl $project_folder/results/MACS2/H3K27ac_IP_merged_peaks.narrowPeak hg38 \
 > $project_folder/results/MACS2/H3K27ac_IP_merged_peaks_annotation.csv \
 -annStats $project_folder/results/MACS2/H3K27ac_IP_merged_peaks_annotation_summary.csv

annotatePeaks.pl $project_folder/results/MACS2/H3K27ac_negIP_merged_peaks.narrowPeak hg38 \
 > $project_folder/results/MACS2/H3K27ac_negIP_merged_peaks_annotation.csv \
 -annStats $project_folder/results/MACS2/H3K27ac_negIP_merged_peaks_annotation_summary.csv

echo  "$BLUE--------- Read quantification ---------$RESET"

featureCounts -T 9 -s 1 -t gene -g gene_id -a $annotation \
 -o $project_folder/results/featureCount/readCount_gene_gencode.v27_SE_wInput.txt $project_folder/results/aligned/*.filtered.bam

echo  "$BLUE--------- Enrichment analysis ---------$RESET"

Rscript -e "rmarkdown::render('$project_folder/scripts/02_ChIRP-seq_enrichment_analysis.Rmd',params=list(args = myarg))"