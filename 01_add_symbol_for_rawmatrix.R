#######change gene ID and get matrix ######
rm(list = ls())
setwd("D:/project/Maternal_gene")
counts<-read.table("./count_matrix/join.count",header=T)
head(counts)
tail(counts)
#rownames(counts) <- counts
library('biomaRt')
library("curl")
library(ggrepel)
library(dplyr)
my_mart <-useMart("ensembl")
list<-listDatasets(my_mart)
#####查看database中有哪些物种的基因####
mart <- useDataset("rnorvegicus_gene_ensembl", 
                   useMart("ensembl"))#####Rat gene:rnorvegicus_gene_ensembl
listAttributes(mart)

my_ensembl_gene_id<-counts$gene
head(my_ensembl_gene_id)
options(timeout = 4000000)
#####use ensembl_transcript_id to trans######
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',
                                 'ensembl_transcript_id',"description"),
                    filters = 'ensembl_gene_id', 
                    values= my_ensembl_gene_id, mart = mart)
mms_symbols$gene<-mms_symbols$ensembl_gene_id
counts<-merge(counts,mms_symbols,by="gene")
str(counts)
counts<-counts[,c(2:18)]
counts<-counts[!duplicated(counts$external_gene_name), ]
rownames(counts)<-counts$external_gene_name
head(counts)
counts<-counts[,c(4:9,1:3,10:15)]
colnames(counts)<-c("AdVector.1","AdVector.2","AdVector.3",
"AdVector.4","AdVector.5","AdVector.6",
"AdG1.1"  ,"AdG1.2"   ,"AdG1.3"  , 
"AdZ4.1"  ,"AdZ4.2"   ,"AdZ4.3"  , 
"AdZ4G1.1","AdZ4G1.2" ,"AdZ4G1.3")

write.table(counts,"./count_matrix/rawcounts_add_genesymbol.txt")
write.csv(counts,"./count_matrix/rawcounts_add_genesymbol.csv")
# global raw_counts pca
library(ggpubr)
library(ggthemes)
library(gmodels)
library(export)
pca.info <- fast.prcomp(counts)
head(pca.info)
summary(pca.info) 
head(pca.info$rotation) 

pca.data <- data.frame(sample = rownames(pca.info$rotation), 
                       Type=c(rep("AdVector_b1",3),rep("AdVector_b2",3),
                              rep("AdG1",3),rep("AdZ4",3),
                              rep("AdZ4G1",3)), pca.info$rotation)
p<-ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=TRUE, ellipse.type="convex")
p

