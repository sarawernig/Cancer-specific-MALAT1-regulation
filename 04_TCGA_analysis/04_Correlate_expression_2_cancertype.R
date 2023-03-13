######### try to download gene expression data ##########
setwd("~/Analysis/06_201903_Triplex_Integrated/NR4A1_peak/Cancer-ATAC-ATLAS/peak_matrix/expression_data/")
library(TCGAutils)
library(TCGAbiolinks)
library(SummarizedExperiment)

## get expression data FPKM-Q75
exp<-readRDS(file="expression_matching_ATAC_rmDup.rds")

################################ get gene #############
wx.nr4a1<-which(rowData(exp)$external_gene_name=="NR4A1")  
ens.nr4a1<-rowData(exp)[wx.nr4a1,"ensembl_gene_id"]
exp.nr4a1<-assay(exp[ens.nr4a1,])
colnames(exp.nr4a1)<-TCGAbarcode(colnames(exp.nr4a1), sample = T)



#### get ATAC peak signal ###
id32<-read.delim("../PRAD_72732_signal_dedup.txt", check.names = F)
id33<-read.delim("../PRAD_72733_signal_dedup.txt", check.names = F)


#### get meta data and cancer types 
meta<-read.delim("../clinical.cases/clinical.tsv")
cancer.types<-split(paste0(meta$submitter_id,"-01A"),meta$project_id)
dir.create("cancer.types_Corr")

for(i in names(cancer.types)){
  print(i)
  id<-sapply(strsplit(i,"-"), function(x) x[2])
  # get barcodes
  select.id<-intersect(cancer.types[[i]],colnames(exp.nr4a1))
  
  # check if data is available
  if(length(select.id)>0) {
    # get data
    df<-data.frame(
      x=log2(as.vector(exp.nr4a1[,select.id])),
      y.72732=as.numeric(id32[select.id]),
      y.72733=as.numeric(id33[select.id])
    )
    
    
    fit.72732<-lm(y.72732~x, df)
    sum.fit.72732<-summary(fit.72732)
    
    fit.72733<-lm(y.72733~x, df)
    sum.fit.72733<-summary(fit.72733)
    pdf(paste0("cancer.types_Corr/",id,"_correlation_NR4A1_72732.pdf"), height = 4, width=5)
    with(df,plot(x,y.72732,pch=19,
                 main=paste0(id," - NR4A1 - peak PRAD_72732"),cex=0.5,
                 xlab="RNA expression (log2(FPKMs))", ylab="ATAC signal log2(CPM)"))
    abline(fit.72732)
    with(df,mtext(paste0("r=",round(sqrt(abs(sum.fit.72732$adj.r.squared)),digits = 2), " - ", "pVal=",sum.fit.72732$coefficients[2,4]),
                  side = 3))
    dev.off()
    
    
    
    pdf(paste0("cancer.types_Corr/",id,"_correlation_NR4A1_72733.pdf"), height = 4, width=5)
    with(df,plot(x,y.72733,pch=19,
                 main=paste0(id," - NR4A1 - peak PRAD_72733"),cex=0.5,
                 xlab="RNA expression (log2(FPKMs))", ylab="ATAC signal log2(CPM)"))
    abline(fit.72733)
    with(df,mtext(paste0("r=",round(sqrt(abs(sum.fit.72733$adj.r.squared)),digits = 2), " - ", "pVal=",sum.fit.72733$coefficients[2,4]),
                  side = 3))
    dev.off()
  }
  
  
  
}



### correlate with ATAC peak #########
id32<-read.delim("PRAD_72732_signal_dedup.txt", check.names = F)
id33<-read.delim("PRAD_72733_signal_dedup.txt", check.names = F)
RNA.codes<-TCGAbarcode(colnames(exp.nr4a1), sample = T)



fit.72732<-lm(y.72732~x, df)
sum.fit.72732<-summary(fit.72732)

fit.72733<-lm(y.72733~x, df)
sum.fit.72733<-summary(fit.72733)


# R-square
sqrt(sum.fit$adj.r.squared)

# r
sqrt(sum.fit$adj.r.squared)

# p-value of slope
sum.fit$coefficients[2,4]




