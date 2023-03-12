---
title: "R Notebook"
output: html_notebook
---

library(rtracklayer)

rtracklayer::import('/01_Projects/genome/gencode.v27.annotation+chrR-Reversed_bin500.gtf') -> GTF
head(GTF)

new_GTF <- as.data.frame(GTF)
new_GTF<- data.frame(new_GTF$gene_id, new_GTF$gene_name,new_GTF$gene_type)

colnames(new_GTF) <- c("GeneID","GeneName","GeneBiotype")

unique(new_GTF$GeneBiotype)

new_GTF$GeneID <- gsub("(ENSG[0-9]+).[0-9]+","\\1",new_GTF$GeneID)

new_GTF <- unique(new_GTF)

head(new_GTF)
nrow(new_GTF)

# Import results table

seq <- read.table("/01_Projects/lncRNA_project/ChIRP-seq/results/Enrichment/ChIRP-seq_results.csv", 
                  sep = ",",header = TRUE)
nrow(seq)
length(unique(annotated$GeneID))

head(seq)

annotated <- merge(seq,new_GTF,by="GeneID",all.x = TRUE)
nrow(annotated)
head(annotated)
res <- as.data.frame(table(annotated$GeneBiotype))
colnames(res) <- c("GeneBiotype","Frequency")
res
nrow(annotated)
length(unique(annotated$GeneID))

write.table(annotated, file="/01_Projects/lncRNA_project/ChIRP-seq/results/Enrichment/ChIRP-seq_results8_anno.full.csv",
            sep="\t",quote = FALSE,row.names = FALSE)

up <- subset(annotated,log2FoldChange > 0)
res.up <- as.data.frame(table(up$GeneBiotype))
colnames(res.up) <- c("GeneBiotype","Frequency")
res.up

write.table(res.up, file="/01_Projects/lncRNA_project/ChIRP-seq/results/Enrichment/ChIRP-seq_results_anno-up.csv",
            sep="\t",quote = FALSE,row.names = FALSE)
