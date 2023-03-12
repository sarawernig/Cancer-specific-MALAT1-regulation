#!/usr/bin/env Rscript

print("DiffBind analysis")

### INSTALLATION ###
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("DiffBind", version = "3.8")

library("DiffBind")

#Make a subdirectory for results
outpath = "/01_Projects/lncRNA_project/ATAC-seq_HeLa_48h/results/diffBind/"
subDir <- "MALAT1vsCTRL"

concat = function(v) {
    resu = ""
    for (i in 1:length(v))
    resu = paste(resu,v[i],sep="")
    resu
}
if (file.exists(subDir)){
    setwd(file.path(outpath, subDir))
} else {
    dir.create(file.path(outpath, subDir))
    setwd(file.path(outpath, subDir))
}
outdir <- concat(c(outpath,subDir,"/"))


samples <- read.csv("/01_Projects/lncRNA_project/ATAC-seq_HeLa_48h/results/diffBind/sampleSheet_diffBind_MALAT1.csv", header=TRUE, sep=";")
samples

DBdata <- dba(sampleSheet=samples)
DBdata

pdf("Correlation_heatmap.pdf")
plot(DBdata)
dev.off()

DBdata <- dba.count(DBdata)
DBdata

DBdata <- dba.contrast(DBdata, group1=DBdata$masks$Normal, group2=DBdata$masks$MALAT1_KD, name1="Control", name2="MALAT1_KD")

dba.show(DBdata, bContrasts=TRUE) #Gives info on the contracts (comparisons)

DBdata <- dba.analyze(DBdata, method=DBA_EDGER)
DBdata <- dba.analyze(DBdata, method=DBA_DESEQ2)
DBdata

pdf("diffBind_heatmap.pdf")
plot(DBdata)
dev.off()

pdf("binding-affinity_heatmap.pdf")
dba.plotHeatmap(DBdata,contrast=1,correlations=FALSE,method=DBA_EDGER)
dev.off()
#
pdf("MAplot.pdf")
dba.plotMA(DBdata,method=DBA_EDGER)
dev.off()

pdf("MAplot_compareMethods.pdf")
par(mfrow=c(3,1))
dba.plotMA(DBdata,method=DBA_EDGER,bNormalized=FALSE)
dba.plotMA(DBdata,method=DBA_EDGER,bNormalized=TRUE)
dba.plotMA(DBdata,method=DBA_DESEQ2,bNormalized=TRUE)
dev.off()

pdf("Volcano_plot.pdf")
dba.plotVolcano(DBdata,method=DBA_EDGER)
dev.off()

pdf("VennDiagram_compareMethods.pdf")
dba.plotVenn(DBdata,contrast=1,method=DBA_ALL_METHODS)
dev.off()

pdf("PCA_plot.pdf")
dba.plotPCA(DBdata,method=DBA_EDGER)
dba.plotPCA(DBdata,method=DBA_EDGER,contrast=1)
dev.off()

pdf("Boxplot.pdf")
dba.plotBox(DBdata,method=DBA_EDGER)
dev.off()

result <- dba.report(DBdata,method=DBA_EDGER)
result

result2 <- dba.report(DBdata,method=DBA_DESEQ2)
result2

write.csv(result,file="Differentially_bound_peaks_edgeR.csv")
write.csv(result2,file="Differentially_bound_peaks_deseq2.csv")

print("The End")

#savehistory(file="diffBind_commands.R")
