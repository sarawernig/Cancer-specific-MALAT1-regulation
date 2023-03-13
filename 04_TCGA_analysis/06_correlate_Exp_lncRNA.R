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

wx.malat<-which(rowData(exp)$external_gene_name=="MALAT1")  
ens.malat<-rowData(exp)[wx.malat,"ensembl_gene_id"]

wx.ha<-which(rowData(exp)$external_gene_name=="HOTAIR")  
ens.ha<-rowData(exp)[wx.ha,"ensembl_gene_id"]

wx.neat1<-which(rowData(exp)$external_gene_name=="NEAT1")  
ens.neat1<-rowData(exp)[wx.neat1,"ensembl_gene_id"]

exp.lncRNA<-assay(exp[c(ens.nr4a1,ens.malat,ens.ha,ens.neat1),])
colnames(exp.lncRNA)<-TCGAbarcode(colnames(exp.lncRNA), sample = T)
rownames(exp.lncRNA)<-c("NR4A1","MALAT1","HOTAIR", "NEAT1")
############################################################################
## get peak data ###
id32<-read.delim("../PRAD_72732_signal_dedup.txt", check.names = F)
id33<-read.delim("../PRAD_72733_signal_dedup.txt", check.names = F)


df<-data.frame(
  x=log2(exp.lncRNA["NR4A1",]+1),
  x.72732=as.numeric(id32[colnames(exp.lncRNA)]),
  x.72733=as.numeric(id33[colnames(exp.lncRNA)]),
  y.malat=log2(exp.lncRNA["MALAT1",]+1),
  y.hotair=log2(exp.lncRNA["HOTAIR",]+1),
  y.neat1=log2(exp.lncRNA["NEAT1",]+1)
)

colnames(df)<-c("NR4A1", "ID_72732","ID_72733","MALAT1","HOTAIR", "NEAT1")
rownames(df)<-colnames(exp.lncRNA)


dir.create("lncRNA_Corr")
path<-"lncRNA_Corr/all"
dir.create(path)


###  plot

for(x in c("NR4A1", "ID_72733", "ID_72732")){
  for(y in c("NEAT1","HOTAIR", "MALAT1")){
    # filter zero counts
    rmv<-rownames(df)[which(df[,y]!=0)]
    fit<-lm(df[rmv,y]~df[rmv,x])
    sum.fit<-summary(fit)
    
    pdf(paste0(path,"/",x,"_",y,"_correlation.pdf"), height = 4, width=5)
    plot(df[rmv,x],df[rmv,y],pch=19,
         main=paste0(x," - ",y),cex=0.5,
         xlab=x, ylab=y)
    abline(fit)
    with(df,mtext(paste0("r=",round(sqrt(abs(sum.fit$adj.r.squared)),digits = 2), " - ", "pVal=",sum.fit$coefficients[2,4]),
                  side = 3))
    dev.off()
  }

}


### go cancer type specific
meta<-read.delim("../clinical.cases/clinical.tsv")
cancer.types<-split(paste0(meta$submitter_id,"-01A"),meta$project_id)

for(i in names(cancer.types)){
  print(i)
  id<-sapply(strsplit(i,"-"), function(x) x[2])
  # get barcodes
  select.id<-intersect(cancer.types[[i]],rownames(df))
  
  if(length(select.id)>0) {
    path<-paste0("lncRNA_Corr/",id)
    dir.create(path)
    
    for(x in c("NR4A1", "ID_72733", "ID_72732")){
      for(y in c("NEAT1","HOTAIR", "MALAT1")){
        # remove zero counts
        select.id.rmv<-select.id[which(df[select.id,y]!=0)]
        if(length(select.id.rmv)>1){
          fit<-lm(df[select.id.rmv,y]~df[select.id.rmv,x])
          sum.fit<-summary(fit)
          
          pdf(paste0(path,"/",x,"_",y,"_",id,"_correlation.pdf"), height = 4, width=5)
          plot(df[select.id.rmv,x],df[select.id.rmv,y],pch=19,
               main=paste0(id," - ",x," - ",y),cex=0.5,
               xlab=x, ylab=y)
          abline(fit)
          with(df,mtext(paste0("r=",round(sqrt(abs(sum.fit$adj.r.squared)),digits = 2), " - ", "pVal=",sum.fit$coefficients[2,4]),
                        side = 3))
          dev.off()
        }
      }
    }
  }

}

