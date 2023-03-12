#!/usr/bin/env Rscript

library("dplyr")
library("ggplot2")
library("ggpubr")
library("tidyverse")
library("reshape2")
library("edgeR")
library("ggpubr")
library("ggpmisc")
library("cowplot")
library("pheatmap")
library("ComplexHeatmap")


setwd("~/01_Projects/lncRNA_project/out/")
load("~/01_Projects/lncRNA_project/out/.RData")

######################## Scatter plot Figure 5, A #############################

atac <- read.delim2(file="panc_ATAC-seq_rawCounts.tab")
rna <- read.delim2(file="GSE124231_counts_rnaseq.txt",row.names = 1,header=T)

counts <- DGEList(counts=as.matrix(rna))
counts <- calcNormFactors(counts)
norm.counts <- cpm(counts)
norm.counts <- as.data.frame(norm.counts)


NR4A1 <- norm.counts %>% filter(row.names(norm.counts) == "ENSG00000123358")


peak1 <- atac[1,]
peak1 <- peak1[,4:ncol(peak1)]

peak2 <- atac[2,]
peak2 <- peak2[,4:ncol(peak2)]

colnames(peak1) <-
  str_replace_all(colnames(peak1), regex("GSM[0-9]{7}_Sample_ATAC_SD_"), "")
colnames(peak1) <-
  str_replace_all(colnames(peak1), regex("_[0-9]{8}"), "")
colnames(peak1)


colnames(peak2) <-
  str_replace_all(colnames(peak2), regex("GSM[0-9]{7}_Sample_ATAC_SD_"), "")
colnames(peak2) <-
  str_replace_all(colnames(peak2), regex("_[0-9]{8}"), "")
colnames(peak2)

peak1[match(row.names(peak1), row.names(NR4A1)),"c"]

peak1 <- t(peak1)
peak1 <- as.data.frame(peak1)


peak2 <- t(peak2)
peak2 <- as.data.frame(peak2)

NR4A1 <- t(NR4A1)
NR4A1 <- as.data.frame(NR4A1)

p1 <- merge(peak1, NR4A1, by=0)
p2 <- merge(peak2, NR4A1, by=0)

colnames(p1) <- c("ID","peak1","NR4A1")
colnames(p2) <- c("ID","peak2","NR4A1")


p1$NR4A1 <- log2(p1$NR4A1)
p2$NR4A1 <- log2(p2$NR4A1)


p1$peak1 <- as.numeric(p1$peak1)
p2$peak2 <- as.numeric(p2$peak2)


p1.lm <- lm(peak1 ~ NR4A1, p1)
p2.lm <- lm(peak2 ~ NR4A1, p2)


plot1 <- ggplot(data = p1, mapping = aes(x = NR4A1,y = peak1))+
  geom_point(size = 3, colour = "#F7766C") +
  stat_smooth(method = 'lm', colour = "black") + theme_cowplot(12) +
  # adds R^2 and p-value
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           r.accuracy = 0.00001,
           p.accuracy = 0.00001,
           size = 4) 

plot2 <- ggplot(data = p2, mapping = aes(x = NR4A1,y = peak2))+
  geom_point(size = 3, colour = "#07BEC4") +
  stat_smooth(method = 'lm', colour = "black") + theme_cowplot(12) +
  # adds R^2 and p-value
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           r.accuracy = 0.00001,
           p.accuracy = 0.00001,
           size = 4) 

#peak1: #F7766C
#peak2: #07BEC4


pdf("PANC1_ROIvsNR4A1.pdf")
plot_grid(plot1, plot2)
dev.off()


svg("PANC1_ROIvsNR4A1.svg")
plot_grid(plot1, plot2)
dev.off()

save.image("~/01_Projects/lncRNA_project/out/.RData")

