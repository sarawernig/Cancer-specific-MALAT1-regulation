setwd("~/Analysis/06_201903_Triplex_Integrated/NR4A1_peak/Cancer-ATAC-ATLAS/peak_matrix/")
library(TCGAutils)
#get meta data
meta<-read.delim("clinical.cases/clinical.tsv")



id33<-read.delim("PRAD_72733_signal_dedup.txt", check.names = F)
id32<-read.delim("PRAD_72732_signal_dedup.txt", check.names = F)
barcodes<-TCGAbarcode(colnames(id32)[c(-1)])


#get.both in graph
df.both<-data.frame("Peak-ID"=rep(c("PRAD_72732","PRAD_72733"), each=length(mx)),
              "ATAC-signal"=c(as.numeric(id32[1,c(-1)]),as.numeric(id33[1,c(-1)])),
               "cancer_type"=rep(gsub("TCGA-","",as.character(meta$project_id[mx])),length(mx)))

# sort by ATAC signal
mean.val<-sapply(split(df.both$ATAC.signal,df.both$cancer_type ),mean)
ord<-order(mean.val)
df.both$cancer_type<-factor(df.both$cancer_type, levels=levels(df.both$cancer_type)[ord])

library(ggplot2)
library(ggbeeswarm)

g<-ggplot(df.both,aes(cancer_type,ATAC.signal, fill=Peak.ID))+
  geom_boxplot()+theme_bw()

pdf("cancerType_distribution.pdf", width=11, height = 4.5)
  print(g)
dev.off()




################################## get single ones ###################################

mx<-match(barcodes,meta$submitter_id)
df<-data.frame("PRAD_72732"=as.numeric(id32[1,c(-1)]),
               "PRAD_72733"=as.numeric(id33[1,c(-1)]),
               "cancer_type"=gsub("TCGA-","",as.character(meta$project_id[mx])))
row.names(df)<-barcodes


df$cancer_type<-factor(df$cancer_type, levels=levels(df.both$cancer_type))


g<-ggplot(df,aes(cancer_type,PRAD_72733))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(bandwidth = 0.25)+theme_bw()

pdf("cancerType_distribution_PRAD_72733.pdf", width=10, height = 4.5)
  print(g)
dev.off()

#######

g<-ggplot(df,aes(cancer_type,PRAD_72732))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(bandwidth = 0.25)+theme_bw()

pdf("cancerType_distribution_PRAD_72732.pdf", width=10, height = 4.5)
  print(g)
dev.off()


