#!/bin/bash

#Run as: bash ./01_RNA-seq_48h_analysis.sh

echo "RNA-seq 48h post-MALAT1-Kd analysis"

project_folder="/01_Projects/lncRNA_project/RNA-seq_HeLa_48h"
genome="/01_Projects/genome/STAR-index-hg38+chrR_reversed-bigOnly"
annotation="/01_Projects/genome/gencode.v27.annotation+chrR-Reversed_bin500.gtf"
rDNA="/01_Projects/genome/rDNA_genBank_U13369.1.fasta"
chromSizes="/01_Projects/genome/hg38_chrR_U13369_bigChr.chrom.sizes"
blacklistedRegions="/01_Projects/genome/GRCh38.blacklist.bed"

cd $project_folder
mkdir raw_data results trimmed logs scripts aligned

cd raw_data7
echo "----- Demultiplexing -----"
bcl2fastq --runfolder-dir . --output-dir fastQ/ --no-lane-splitting --sample-sheet sampleSheet.txt --minimum-trimmed-read-length 35 --mask-short-adapter-reads 0 -p 12 --use-bases-mask Y*,I8Y*,Y*

echo "----- FastQC analysis -----"
mkdir $project_folder/results/fastQC/
fastqc $project_folder/raw_data/*.fastq.gz -o $project_folder/results/fastQC/ -t 10
multiqc $project_folder/results/fastQC/ -o $project_folder/results/fastQC/

echo "----- Ovation SOLO specific trimming (first 5 bases of forward reads) -----"
for forward in $project_folder/raw_data/*_R1_001.fastq.gz; do cutadapt -u 5 -o ${forward/%.fastq.gz/.trim.fastq.gz} $forward; done

echo "----- NextSeq550 specific trimming of polyG/polyA/polyT tails from reverse reads -----"
for reverse in $project_folder/raw_data/*_R2_001.fastq.gz; do cutadapt -a "A{100}" -a "G{100}" -a "T{100}" -o ${reverse/%.fastq.gz/.trim.fastq.gz} $reverse; done

echo "----- Adapter trimming, quality filtering and read pairing -----"
for forward in $project_folder/raw_data/*_R1_001.trim.fastq.gz ; do java -jar trimmomatic-0.38.jar PE \
-phred33 -threads 20 $forward ${forward/%_R1_001.trim.fastq.gz/_R2_001.trim.fastq.gz} ${forward/%_R1_001.trim.fastq.gz/_1_paired.fastq} \
${forward/%_R1_001.trim.fastq.gz/_1_unpaired.fastq} ${forward/%_R1_001.trim.fastq.gz/_2_paired.fastq} ${forward/%_R1_001.trim.fastq.gz/_2_unpaired.fastq} \
ILLUMINACLIP:/Users/sara/tools/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:30 SLIDINGWINDOW:4:20; done

echo "----- FastQC analysis after trimming and pairing -----"
fastqc $project_folder/raw_data/*_paired* -t 20 -o $project_folder/results/fastQC/paired/
multiqc $project_folder/results/fastQC/paired/ -o $project_folder/results/fastQC/paired/

echo "----- Remove ribosomal RNA -----"
for file in $project_folder/raw_data/*_1_paired.fastq; do bbduk.sh -Xmx24g minlength=30 \
in=$file in2=${file/%_1_paired.fastq/_2_paired.fastq} out=${file/%_1_paired.fastq/_rRNA_removed_1.fastq} \
out2=${file/%_1_paired.fastq/_rRNA_removed_2.fastq} ref=$rDNA stats=logs/${file/%_1_paired.fastq/_stats.txt}; done

echo "----- FastQC analysis after rRNA removal -----"
mkdir $project_folder/results/fastQC/rRNA_removed
fastqc $project_folder/raw_data/*_rRNA_removed* -t 20 -o $project_folder/results/fastQC/rRNA_removed
multiqc $project_folder/results/fastQC/rRNA_removed -o $project_folder/results/fastQC/rRNA_removed

echo "----- Alignment to the reference genome -----"
cd $project_folder/raw_data
for file in *_rRNA_removed_1.fastq; do STAR --runMode alignReads \
--runThreadN 15 --sjdbGTFfile $annotation \
--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1  \
--outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--outSAMtype BAM SortedByCoordinate --outWigType bedGraph --genomeDir /Users/sara/data/hg38_genome/STAR_index_38bp/ \
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 \
-readFilesIn $file ${file/%_1.fastq/_2.fastq} --outFileNamePrefix ${file/%_rRNA_removed_1.fastq/_}; done

echo "----- Quality filtering -----"
for file in $project_folder/aligned/*_Aligned.sortedByCoord.out.bam; do samtools view -@ 20 -b -F 4 -q 20 $file > \
${file/%_Aligned.sortedByCoord.out.bam/_filtered.bam}; done

echo "----- Mapping quality check -----"
cd ..
mkdir $project_folder/results/qualimap
mv $project_folder/raw_data/*bam aligned/
cd $project_folder/aligned
for file in *_filtered.bam; do qualimap rnaseq -p strand-specific-forward -gtf $annotation --java-mem-size=16G \
-pe -bam $file -outformat PDF -outdir ../results/qualimap/${file/%_filtered.bam/}; done

echo "----- Convert RNAseq BAM files to normalized TDF files -----"

for file in *_filtered.bam; do samtools1.9 index -b $file; done
for file in *_filtered.bam; do bamCoverage -b $file -o ${file/%.bam/_reverseStrand.bigwig} --filterRNAstrand forward \
--binSize 25 --blackListFileName $blacklistedRegions -p 20 --normalizeUsing RPKM --minMappingQuality 30; done

for file in *_filtered.bam; do bamCoverage -b $file -o ${file/%.bam/_forwardStrand.bigwig} --filterRNAstrand reverse \
--binSize 25 --blackListFileName $blacklistedRegions -p 20 --normalizeUsing RPKM --minMappingQuality 30; done

for file in aligned/*.bigwig; do ~/tools/UCSC_tools/bigWigToWig $file ${file/.bigwig/.wig};done

for file in aligned/*.wig; do igvtools toTDF $file ${file/%.wig/.tdf} $chromSizes; done

echo "----- READ COUNT -----"
cd ..
mkdir $project_folder/results/featureCount
featureCounts -p -T 20 -s 1 -t gene -g gene_id -C -B -a $annotation -o $project_folder/results/featureCount/readCount_gene_PE.txt $project_folder/aligned/*_filtered.bam

echo "----- DIFFERENATIAL EXPRESSION ANALYSIS -----"
mkdir $project_folder/results/DEG_analysis

Rscript -e "rmarkdown::render('$project_folder/scripts/02_RNA-seq_48h_DEGs.Rmd',params=list(args = myarg))"





