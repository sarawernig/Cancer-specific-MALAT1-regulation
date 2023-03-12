#!/bin/bash

set -e

echo "----- ATAC-seq PIPELINE (PE only) -----"
annotation="/01_Projects/genome/gencode.v27.annotation+chrR-Reversed_bin500.gtf"
project_folder="/01_Projects/lncRNA_project/ATAC-seq_HeLa_48h"
genome="/01_Projects/genome/Bowtie2_hg38+chrR-rev_index/hg38+chrR"
fasta="/01_Projects/genome/Hg38_chrR-Reversed_bigOnly.fa"
blacklistedRegions="/01_Projects/genome/GRCh38.blacklist.bed"
chromSizes="/01_Projects/genome/hg38_chrR_U13369_bigChr.chrom.sizes"

cd $project_folder

samples=$project_folder/raw_data/*.fastq.gz

echo "The project folder: " $project_folder
echo "Annotation file: " $annotation
echo "Genome: " $genome
echo "Genome fasta file: " $fasta
echo "Samples: " $(basename ${samples/%.fastq.gz/})

mkdir -p raw_data aligned results logs scripts trimmed

echo "----- Quality Control -----"
mkdir -p $project_folder/results/fastQC/raw/
fastqc --contaminants /Users/sara/tools/Trimmomatic-0.38/adapters/NGS_contaminants.fa --threads 20 --outdir $project_folder/results/fastQC/raw/ $project_folder/raw_data/*.fastq.gz

multiqc $project_folder/results/fastQC/raw/ -o $project_folder/results/fastQC/raw/

echo "----- Adapter trimming and read pairing -----"
for forward in raw_data/*_R1.fastq.gz; do java -jar /Users/sara/tools/Trimmomatic-0.38/trimmomatic-0.38.jar PE \
-phred33 -threads 18 $forward ${forward/%_R1.fastq.gz/_R2.fastq.gz} ${forward/%_R1.fastq.gz/_R1_paired.fastq} \
 ${forward/%_R1.fastq.gz/_R1_unpaired.fastq} ${forward/%_R1.fastq.gz/_R2_paired.fastq}  ${forward/%_R1.fastq.gz/_R2_unpaired.fastq} \
 ILLUMINACLIP:/Users/sara/tools/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30 2>> logs/Trimmomatic_stats.txt; done

echo "----- Quality control after adapter removal -----"
mkdir -p $project_folder/results/fastQC/trimmed/
mv $project_folder/raw_data/*paired*  $project_folder/trimmed/

fastqc --threads 20 --outdir $project_folder/results/fastQC/trimmed/ $project_folder/trimmed/*_paired*
multiqc $project_folder/results/fastQC/trimmed/ -o $project_folder/results/fastQC/trimmed/

echo "----- Mapping to the reference genome -----"
for file in $project_folder/trimmed/*_R1_paired.fastq; do bowtie2 --threads 16 --local -t --very-sensitive-local -X 5000 \
 --no-discordant --dovetail -x $genome -1 $file -2 ${file/%_R1_paired.fastq/_R2_paired.fastq}  | samtools view -bS -@ 18 -h - > ${file/%_R1_paired.fastq/_PE_aligned.bam}; done

mv $project_folder/trimmed/*.bam $project_folder/aligned/

for file in $project_folder/aligned/*.bam; do samtools sort -@ 12 $file -o ${file/%.bam/.sorted.bam}; done
for file in $project_folder/aligned/*.sorted.bam; do samtools index -b $file; done

echo "----- Mapping quality check -----"
mkdir -p results/qualimap
#how to make the sample_info.txt
for file in $project_folder/trimmed/*1_paired*; do echo ${file/%_R1_paired.fastq/_PE_aligned.bam}; done

qualimap multi-bamqc --java-mem-size=16G -outdir $project_folder/results/qualimap -d $project_folder/aligned/sample_info.txt -r -gff $annotation -outformat PDF

echo "----- Quality filtering -----"
for file in aligned/*.sorted.bam; do samtools view -@ 20 -b -F 4 -q 30 $file > ${file/%.bam/_filtered.bam}; done

echo "----- Remove chrM -----"
for file in aligned/*_aligned.sorted_filtered.bam; do samtools index $file; done

for file in aligned/*_aligned.sorted_filtered.bam; do samtools idxstats -@ 6 $file | cut -f 1 | grep -v chrM | xargs samtools view -@ 16 -b $file \
 > ${file/%_aligned.sorted_filtered.bam/_aligned_filtered_rmChrM.bam}; done

echo "----- Remove duplicates -----"
for file in aligned/*_rmChrM.bam; do java -jar -Xmx32g /Users/sara/tools/picard.jar MarkDuplicates INPUT=$file OUTPUT=${file/%.bam/_noDup.bam} \
 METRICS_FILE=${file/%.bam/_noDup_metrics.txt} REMOVE_DUPLICATES=true CREATE_INDEX=true ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
 USE_JDK_DEFLATER=true USE_JDK_INFLATER=true CREATE_MD5_FILE=true; done


echo "----- Shift reads based on TN5 cuts -----"
for file in aligned/*/*_noDup.bam; do echo $file; alignmentSieve -b $file -o ${file/%.bam/.ATACshift.bam} -p 14 \
 --ATACshift --blackListFileName $blacklistedRegions ; done

echo "----- Make bigwig files for IGV -----"
for file in aligned/*_PE_aligned_filtered_rmChrM_noDup.ATACshift.bam; do echo $file; samtools sort -@ 20 $file > ${file/%.bam/.sort.bam}; done

for file in aligned/*_PE_aligned_filtered_rmChrM_noDup.ATACshift.sort.bam; do echo $file; samtools index -b $file; done

for file in aligned/*_PE_aligned_filtered_rmChrM_noDup.ATACshift.sort.bam; do bamCoverage -b $file -o ${file/%.bam/.bigwig} -p 5 \
 --effectiveGenomeSize 2747877777 --normalizeUsing RPGC; done

for file in aligned/*.ATACshift.sort.bigwig; do echo $file;  ~/tools/UCSC_tools/bigWigToWig $file ${file/.bigwig/.wig}; done

for file in aligned/*.wig; do igvtools toTDF $file ${file/%.wig/.tdf} $chromSizes; done

echo "----- Insert size plot -----"
bamPEFragmentSize --histogram fragmentSize_histogram.png --plotTitle "Fragment size (bp) for HeLa ATAC-seq" --maxFragmentLength 1000 \
 --bamfiles aligned/*_PE_aligned_filtered_rmChrM_noDup.bam --samplesLabel CTRL_rep1 CTRL_rep2 CTRL_rep3 MALAT1_rep1 MALAT1_rep2 MALAT1_rep3

echo "----- Peak calling -----"
mkdir -p results/MACS2
for file in aligned/*_PE_aligned_filtered_rmChrM_noDup.bam; do MACS2 callpeak -f BAMPE -t $file -g hs \
 -n $(basename ${file/%_PE_aligned_filtered_rmChrM_noDup.bam/}) -q 0.1 --shift -100 --extsize 200 -B --SPMR \
 --nomodel --outdir results/MACS2 --keep-dup all --buffer-size=100; done

for file in $project_folder/results/MACS2/*.narrowPeak; do mv $file ${file/%.narrowPeak/.bed}; done

echo "----- Differential accessability after KD: diffBind analysis -----"
mkdir -p $project_folder/results/diffBind
$project_folder/scripts/02_ATAC-seq_48h_diffBind_script_MALAT1.R


echo "----- Peak calling for open regions -----"
mkdir -p $project_folder/results/MACS2/open_regions
for file in $project_folder/aligned/*_PE_aligned_filtered_rmChrM_noDup.bam; do MACS2 callpeak -f BAMPE -t $file -g hs \
 -n $(basename ${file/%_PE_aligned_filtered_rmChrM_noDup.bam/_broad}) --broad --outdir $project_folder/results/MACS2/open_regions \
 --keep-dup all --buffer-size=100; done

for file in $project_folder/results/MACS2/open_regions/*.broadPeak; do mv $file ${file/%.broadPeak/.bed}; done

echo "----- PCA analysis -----"
mkdir -p $project_folder/results/PCA_analysis/
cd $project_folder/aligned
multiBamSummary bins -p 20 --bamfiles Ctrl_rep1_PE_aligned_filtered_rmChrM_noDup.bam Ctrl_rep2_PE_aligned_filtered_rmChrM_noDup.bam Ctrl_rep3_PE_aligned_filtered_rmChrM_noDup.bam \
 MALAT1_rep1_PE_aligned_filtered_rmChrM_noDup.bam MALAT1_rep2_PE_aligned_filtered_rmChrM_noDup.bam MALAT1_rep3_PE_aligned_filtered_rmChrM_noDup.bam \
 -o $project_folder/results/PCA_analysis/multiBamSummary_CTRLvsMALAT1.npz

cd $project_folder/
plotPCA -in results/PCA_analysis/multiBamSummary_CTRLvsMALAT1.npz \
 -o results/PCA_analysis/PCA_readCounts_CTRLvsMALAT1.png -T "PCA of 48h-post MALAT1 KD ATAC-seq samples" --labels Ctrl_rep1 Ctrl_rep2 Ctrl_rep3 MALAT1_rep1 MALAT1_rep2 MALAT1_rep3




