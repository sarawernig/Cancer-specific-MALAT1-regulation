setwd("~/Analysis/06_201903_Triplex_Integrated/NR4A1_peak/Cancer-ATAC-ATLAS/bigwigs/")
library(TCGAbiolinks)

########## download bigwigs #######################################
query <- GDCquery_ATAC_seq(file.type = "bigWigs",tumor = "THCA") 
GDCdownload(query,method = "client")

is(query$results)


######### try to download gene expression data ##########
setwd("~/Analysis/06_201903_Triplex_Integrated/NR4A1_peak/Cancer-ATAC-ATLAS/peak_matrix/")
library(TCGAutils)
id32<-read.delim("PRAD_72732_signal_dedup.txt", check.names = F)
barcodes<-TCGAbarcode(colnames(id32)[c(-1)])

meta<-read.delim("clinical.cases/clinical.tsv")
project.ids<-unique(as.character(meta$project_id))


query <- GDCquery(project = project.ids,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type =  "HTSeq - FPKM-UQ",
                  experimental.strategy = "RNA-Seq",
                  barcode = barcodes,
                  legacy = F)
GDCdownload(query, method = "api", files.per.chunk = 10)
data <- GDCprepare(query)

# get matching samples
codes<-TCGAbarcode(colnames(data), sample = T)
data.match<-data[,which(codes %in% colnames(id32)[c(-1)])]

saveRDS(data.match, file="expression_data/expression_matching_ATAC.rds")

#remove duplitcated barcodes 
barCode.short<-TCGAbarcode(colnames(data.match), sample = T)
wx.dup<-which(barCode.short %in% barCode.short[(duplicated(barCode.short))])
data.match2<-data.match[,-wx.dup]

saveRDS(data.match2, file="expression_data/expression_matching_ATAC_rmDup.rds")

################################ get gene #############
wx.nr4a1<-which(rowData(data.match2)$external_gene_name=="NR4A1")  
ens.nr4a1<-rowData(data.match2)[wx.nr4a1,"ensembl_gene_id"]
exp.nr4a1<-assay(data.match2[ens.nr4a1,])


### correlate with ATAC peak #########
id32<-read.delim("PRAD_72732_signal_dedup.txt", check.names = F)
id33<-read.delim("PRAD_72733_signal_dedup.txt", check.names = F)
RNA.codes<-TCGAbarcode(colnames(exp.nr4a1), sample = T)

df<-data.frame(
  x=log2(as.vector(exp.nr4a1)),
  y.72732=as.numeric(id32[RNA.codes]),
  y.72733=as.numeric(id33[RNA.codes])
)

fit.72732<-lm(y.72732~x, df)
sum.fit.72732<-summary(fit.72732)

fit.72733<-lm(y.72733~x, df)
sum.fit.72733<-summary(fit.72733)


# # R-square
# sqrt(sum.fit$adj.r.squared)
# 
# # r
# sqrt(sum.fit$adj.r.squared)
# 
# # p-value of slope
# sum.fit$coefficients[2,4]



pdf("expression_data/correlation_NR4A1_72732.pdf", height = 4, width=5)
  with(df,plot(x,y.72732,pch=19,
       main="NR4A1 - peak PRAD_72732",cex=0.5,
       xlab="RNA expression (log2(FPKMs))", ylab="ATAC signal log2(CPM)"))
  abline(fit.72732)
  with(df,text(min(x),max(y.72732),
      labels = paste0("r=",round(sqrt(sum.fit.72732$adj.r.squared),digits = 2)), adj=0))
  with(df,text(min(x),max(y.72732)-1,
               labels = paste0("pVal=",sum.fit.72732$coefficients[2,4]), adj=0))
dev.off()




pdf("expression_data/correlation_NR4A1_72733.pdf", height = 4, width=5)
  with(df,plot(x,y.72733,pch=19,
               main="NR4A1 - peak PRAD_72733",cex=0.5,
               xlab="RNA expression (log2(FPKMs))", ylab="ATAC signal log2(CPM)"))
  abline(fit.72733)
  with(df,text(min(x),max(y.72733),
               labels = paste0("r=",round(sqrt(sum.fit.72733$adj.r.squared),digits = 2)), adj=0))
  with(df,text(min(x),max(y.72733)-1,
               labels = paste0("pVal=",sum.fit.72733$coefficients[2,4]), adj=0))
dev.off()

