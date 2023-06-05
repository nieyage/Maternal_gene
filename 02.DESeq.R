library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")
rm(list = ls())
rm(list = ls())
setwd("D:/project/Maternal_gene")
counts<-read.table("./count_matrix/rawcounts_add_genesymbol.txt",header=T)
head(counts)
#remove batch effect#
condition <- factor(c(rep("AdVector",6),rep("AdG1",3),rep("AdZ4",3),rep("AdZ4G1",3)))
batch <- factor(c(rep("batch1",3),rep("batch2",3),rep("batch1",3),rep("batch2",3),rep("batch1",3)))
colData <- data.frame(row.names=colnames(counts),condition=condition,batch=batch)
colData
colData$condition <- relevel(colData$condition, ref="AdVector")
countData <- counts[apply(counts, 1, sum) > 1,] 
countData<-na.omit(countData)
dds<-DESeqDataSetFromMatrix(countData,colData,formula(~condition+batch)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
dds <- DESeq(dds)
res <- results(dds)
#PCA not remove batch
library(ggplot2)
vsd <- vst(dds)
head(vsd)
pca_data <- plotPCA(vsd, intgroup=c("batch","condition"), returnData=T, ntop=500)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p<-ggplot(pca_data, aes(PC1, PC2, shape =batch,color=condition)) +
  geom_point(size=3) +
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
geom_text(aes(label=name),vjust=2)
pdf("./PCA-not_remove_batch.pdf",width = 8,height = 5)
p+theme_bw()
dev.off()
## PCA remove batch 
library(limma)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), c(colData$batch))
pca <- plotPCA(vsd, intgroup=c("batch","condition"), returnData=T, ntop=500)
percentVar <- round(100 * attr(pca, "percentVar"))
p<-ggplot(pca, aes(PC1, PC2, shape =batch,color=condition)) +
  geom_point(size=3) +
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
geom_text(aes(label=name),vjust=2)
pdf("./PCA-remove_batch.pdf",width = 8,height = 5)
p+theme_bw()
dev.off()

 head(assay(vsd))
pcacount<-assay(vsd)
head(pcacount)
write.table(pcacount,"RNAseq-PCAcount.csv")
AdG1_vs_AdVector<-results(dds,contrast = c("condition","AdG1","AdVector"))
AdZ4_vs_AdVector<-results(dds,contrast = c("condition","AdZ4","AdVector"))
AdZ4G1_vs_AdVector<-results(dds,contrast = c("condition","AdZ4G1","AdVector"))

# AdG1_vs_AdVector

AdG1_vs_AdVector_diff_gene <-rownames(subset(AdG1_vs_AdVector,padj < 0.05 & abs(log2FoldChange)> 1))
write.csv(AdG1_vs_AdVector,"AdG1_vs_AdVector.csv")
#volcano plot#
pdf("AdG1_vs_AdVector-volcano.pdf")
data<-as.data.frame(AdG1_vs_AdVector)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'AdG1','AdVector'),
                     'No significant')
data_gene<-as.data.frame(table(data$change))
AdG1_gene<-data_gene[data_gene$Var1=="AdG1",]$Freq
AdVector_gene<-data_gene[data_gene$Var1=="AdVector",]$Freq

p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$padj), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("red", "blue","grey"))+
  xlim(c(-10,10)) +
  ylim(c(0,200))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (padj)",
       title=paste("AdG1_vs_AdVector DEG","\n","AdG1:",AdG1_gene,"genes;", "AdVector:",AdVector_gene,"genes;",sep="")
       )+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

library(pheatmap)
#Heatmap for DEG and GO KEGG #
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%AdG1_vs_AdVector_diff_gene),]
#targetcount<-targetcount[,1:6]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
pdf("AdG1_vs_AdVector_diff_gene-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()


library(clusterProfiler)
library(STRINGdb) 
library(tidyverse) 
library(org.Rn.eg.db) 
library(igraph) 
library(ggraph)
string_db <- STRINGdb$new(version="11.5",species=10116)
#GO and KEGG 
data<-as.data.frame(AdG1_vs_AdVector)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'AdG1','AdVector'),
                     'No significant')
res_data<- data
for (i in c("AdG1","AdVector")){
  print(i);
  gene_in_cluster<- rownames(res_data[res_data$change==i,])
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"),OrgDb = org.Rn.eg.db)
  #GO
  pdf(paste0("AdG1_vs_AdVector--up_in_",i,"_GO.pdf",sep=""))
  for (type in c("BP","CC","MF")){
    print(type);
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = type,
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE);
    data<-as.data.frame(ego)
    if(!nrow(data)==0){
      p<-barplot(ego, showCategory=20);
      print(p);
      write.csv(ego,paste0("AdG1_vs_AdVector--up_in_",i,"_GO_",type,".csv"))
      }
  }
  dev.off()
  #KEGG 
  pdf(paste0("AdG1_vs_AdVector--up_in_",i,"_KEGG.pdf",sep=""))
  ego <- enrichKEGG(
    gene = gene.df$ENTREZID,
    keyType = "kegg",
    organism  = "rno",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05)
  data<-as.data.frame(ego)
  if(!nrow(data)==0){
    p<-barplot(ego, showCategory=20);
    print(p);
    edox<-setReadable(ego,OrgDb=org.Rn.eg.db,keyType = "ENTREZID")
    write.table(edox,paste0("AdG1_vs_AdVector--up_in_",i,"_KEGG_",type,".csv"))
  }
  dev.off()
  # Gene-Concept Network
  pdf(paste0("AdG1_vs_AdVector--up_in_",i,"_Gene-Concept-Network.pdf",sep=""))
  input<-AdG1_vs_AdVector$log2FoldChange
  names(input)<-rownames(AdG1_vs_AdVector)
  p1 <- cnetplot(edox, foldChange=input, cex_label_category=0.8)
  p3 <- cnetplot(edox, foldChange=input, circular = TRUE, colorEdge = TRUE) 
  print(p1);
  print(p3);
  dev.off()
}
#GSEA 
library(enrichplot)
gene=bitr(rownames(AdG1_vs_AdVector),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Rn.eg.db") 

gene_df <- data.frame(logFC=AdG1_vs_AdVector$log2FoldChange, 
                      SYMBOL = rownames(AdG1_vs_AdVector)) 
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL
geneList=sort(geneList,decreasing = T) 
kegmt<-read.gmt("D:/project/HMMR/HMMR-OE-P1/rattus_norvegicus_brown-rat_gmt2.gmt") 
names(geneList)<-toupper(names(geneList))
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.45,
           TERM2GENE = kegmt) 
write.csv(KEGG,"AdG1_vs_AdVector_GSVA.csv")
pdf("AdG1_vs_AdVector_GSVA-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="padj",split=".sign")+facet_wrap(~.sign,scales = "free") 
dev.off()
pdf("AdG1_vs_AdVector_GSVA.pdf",width = 10,height = 10)
for (i in 1:200){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) 
dev.off()
library(STRINGdb) 
library(tidyverse) 
library(clusterProfiler) 
library(org.Rn.eg.db) 
library(igraph) 
library(ggraph)
string_db <- STRINGdb$new( version="11.5",species=10116)   
STRINGdb$methods() 
for (i in c("AdG1","AdVector")){
  print(i);
  gene_in_cluster<- rownames(res_data[res_data$change==i,])
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"),OrgDb = org.Rn.eg.db)
  #GO
  pdf(paste0("AdG1_vs_AdVector--up_in_",i,"-Stringdb-network.pdf",sep=""))
  {data_mapped <- gene.df %>% string_db$map(my_data_frame_id_col_names = "ENTREZID",removeUnmappedRows = TRUE)  
    string_db$plot_network( data_mapped$STRING_id )
    data_links <- data_mapped$STRING_id %>% string_db$get_interactions() 
    links <- data_links %>%
      mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
      mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
      dplyr::select(from, to , last_col()) %>% 
      dplyr::rename(weight = combined_score)
    nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
    net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
    igraph::V(net)$deg <- igraph::degree(net) # 每个节点连接的节点数
    igraph::V(net)$size <- igraph::degree(net)/5 #
    igraph::E(net)$width <- igraph::E(net)$weight/10
    ggraph(net,layout = "kk")+
      geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
      geom_node_point(aes(size=size), color="orange", alpha=0.7)+
      geom_node_text(aes(filter=deg>1, label=name), size = 5, repel = T)+
      scale_edge_width(range = c(0.2,1))+
      scale_size_continuous(range = c(1,10) )+
      guides(size=F)+
      theme_graph()
    dev.off()
  }
}


# AdZ4_vs_AdVector
AdZ4_vs_AdVector_diff_gene <-rownames(subset(AdZ4_vs_AdVector,padj < 0.05 & abs(log2FoldChange)> 1))
write.csv(AdZ4_vs_AdVector,"AdZ4_vs_AdVector.csv")
#volcano plot#
pdf("AdZ4_vs_AdVector-volcano.pdf")
data<-as.data.frame(AdZ4_vs_AdVector)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'AdZ4','AdVector'),
                     'No significant')
data_gene<-as.data.frame(table(data$change))
AdZ4_gene<-data_gene[data_gene$Var1=="AdZ4",]$Freq
AdVector_gene<-data_gene[data_gene$Var1=="AdVector",]$Freq
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$padj), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("red", "blue","grey"))+
  xlim(c(-10,10)) +
  ylim(c(0,200))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (padj)",
       title=paste("AdZ4_vs_AdVector DEG","\n","AdZ4:",AdZ4_gene,"genes;", "AdVector:",AdVector_gene,"genes;",sep="")
       )+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()
library(pheatmap)
#Heatmap for DEG and GO KEGG #
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%AdZ4_vs_AdVector_diff_gene),]
#targetcount<-targetcount[,1:6]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
pdf("AdZ4_vs_AdVector_diff_gene-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
data<-as.data.frame(AdZ4_vs_AdVector)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'AdZ4','AdVector'),
                     'No significant')

res_data<- data 
#GO and KEGG 
for (i in c("AdZ4","AdVector")){
  print(i);
  gene_in_cluster<- rownames(res_data[res_data$change==i,])
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"),OrgDb = org.Rn.eg.db)
  #GO
  pdf(paste0("AdZ4_vs_AdVector--up_in_",i,"_GO.pdf",sep=""))
  for (type in c("BP","CC","MF")){
    print(type);
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = type,
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE);
    data<-as.data.frame(ego)
    if(!nrow(data)==0){
      p<-barplot(ego, showCategory=20);
      print(p);
      write.csv(ego,paste0("AdZ4_vs_AdVector--up_in_",i,"_GO_",type,".csv"))
      }
  }
  dev.off()
  #KEGG 
  pdf(paste0("AdZ4_vs_AdVector--up_in_",i,"_KEGG.pdf",sep=""))
  ego <- enrichKEGG(
    gene = gene.df$ENTREZID,
    keyType = "kegg",
    organism  = "rno",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05)
  data<-as.data.frame(ego)
  if(!nrow(data)==0){
    p<-barplot(ego, showCategory=20);
    print(p);
    edox<-setReadable(ego,OrgDb=org.Rn.eg.db,keyType = "ENTREZID")
    write.table(edox,paste0("AdZ4_vs_AdVector--up_in_",i,"_KEGG_",type,".csv"))
  }
  dev.off()
  # Gene-Concept Network
  pdf(paste0("AdZ4_vs_AdVector--up_in_",i,"_Gene-Concept-Network.pdf",sep=""))
  input<-AdZ4_vs_AdVector$log2FoldChange
  names(input)<-rownames(AdZ4_vs_AdVector)
  p1 <- cnetplot(edox, foldChange=input, cex_label_category=0.8)
  p3 <- cnetplot(edox, foldChange=input, circular = TRUE, colorEdge = TRUE) 
  print(p1);
  print(p3);
  dev.off()
}
#GSEA 
library(enrichplot)
gene=bitr(rownames(AdZ4_vs_AdVector),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Rn.eg.db") 

gene_df <- data.frame(logFC=AdZ4_vs_AdVector$log2FoldChange, 
                      SYMBOL = rownames(AdZ4_vs_AdVector)) 
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL
geneList=sort(geneList,decreasing = T) 
kegmt<-read.gmt("D:/project/HMMR/HMMR-OE-P1/rattus_norvegicus_brown-rat_gmt2.gmt") 
names(geneList)<-toupper(names(geneList))
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.45,
           TERM2GENE = kegmt) 
write.csv(KEGG,"AdZ4_vs_AdVector_GSVA.csv")
pdf("AdZ4_vs_AdVector_GSVA-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="padj",split=".sign")+facet_wrap(~.sign,scales = "free") 
dev.off()
pdf("AdZ4_vs_AdVector_GSVA.pdf",width = 10,height = 10)
for (i in 1:200){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) 
dev.off()
library(STRINGdb) 
library(tidyverse) 
library(clusterProfiler) 
library(org.Rn.eg.db) 
library(igraph) 
library(ggraph)
string_db <- STRINGdb$new( version="11.5",species=10116)   
STRINGdb$methods() 
for (i in c("AdZ4","AdVector")){
  print(i);
  gene_in_cluster<- rownames(res_data[res_data$change==i,])
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"),OrgDb = org.Rn.eg.db)
  #GO
  pdf(paste0("AdZ4_vs_AdVector--up_in_",i,"-Stringdb-network.pdf",sep=""))
  {data_mapped <- gene.df %>% string_db$map(my_data_frame_id_col_names = "ENTREZID",removeUnmappedRows = TRUE)  
    string_db$plot_network( data_mapped$STRING_id )
    data_links <- data_mapped$STRING_id %>% string_db$get_interactions() 
    links <- data_links %>%
      mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
      mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
      dplyr::select(from, to , last_col()) %>% 
      dplyr::rename(weight = combined_score)
    nodes <- links %>% {data.frame(gene = c(.$from, .$to)) } %>% distinct()
    net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
    igraph::V(net)$deg <- igraph::degree(net) # 每个节点连接的节点数
    igraph::V(net)$size <- igraph::degree(net)/5 #
    igraph::E(net)$width <- igraph::E(net)$weight/10
    ggraph(net,layout = "kk")+
      geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
      geom_node_point(aes(size=size), color="orange", alpha=0.7)+
      geom_node_text(aes(filter=deg>1, label=name), size = 5, repel = T)+
      scale_edge_width(range = c(0.2,1))+
      scale_size_continuous(range = c(1,10) )+
      guides(size=F)+
      theme_graph()
    dev.off()
  }
}


# AdZ4G1_vs_AdVector

AdZ4G1_vs_AdVector_diff_gene <-rownames(subset(AdZ4G1_vs_AdVector,padj < 0.05 & abs(log2FoldChange)> 1))
write.csv(AdZ4G1_vs_AdVector,"AdZ4G1_vs_AdVector.csv")
#volcano plot#
pdf("AdZ4G1_vs_AdVector-volcano.pdf")
data<-as.data.frame(AdZ4G1_vs_AdVector)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'AdZ4G1','AdVector'),
                     'No significant')
data_gene<-as.data.frame(table(data$change))
AdZ4G1_gene<-data_gene[data_gene$Var1=="AdZ4G1",]$Freq
AdVector_gene<-data_gene[data_gene$Var1=="AdVector",]$Freq
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$padj), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("red", "blue","grey"))+
  xlim(c(-10,10)) +
  ylim(c(0,200))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (padj)",
       title=paste("AdZ4G1_vs_AdVector DEG","\n","AdZ4G1:",AdZ4G1_gene,"genes;", "AdVector:",AdVector_gene,"genes;",sep="")
       )+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()
library(pheatmap)
#Heatmap for DEG and GO KEGG #
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%AdZ4G1_vs_AdVector_diff_gene),]
#targetcount<-targetcount[,1:6]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
pdf("AdZ4G1_vs_AdVector_diff_gene-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
data<-as.data.frame(AdZ4G1_vs_AdVector)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'AdZ4G1','AdVector'),
                     'No significant')

res_data<- data 
#GO and KEGG 
for (i in c("AdZ4G1","AdVector")){
  print(i);
  gene_in_cluster<- rownames(res_data[res_data$change==i,])
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"),OrgDb = org.Rn.eg.db)
  #GO
  pdf(paste0("AdZ4G1_vs_AdVector--up_in_",i,"_GO.pdf",sep=""))
  for (type in c("BP","CC","MF")){
    print(type);
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = type,
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE);
    data<-as.data.frame(ego)
    if(!nrow(data)==0){
      p<-barplot(ego, showCategory=20);
      print(p);
      write.csv(ego,paste0("AdZ4G1_vs_AdVector--up_in_",i,"_GO_",type,".csv"))
      }
  }
  dev.off()
  #KEGG 
  pdf(paste0("AdZ4G1_vs_AdVector--up_in_",i,"_KEGG.pdf",sep=""))
  ego <- enrichKEGG(
    gene = gene.df$ENTREZID,
    keyType = "kegg",
    organism  = "rno",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05)
  data<-as.data.frame(ego)
  if(!nrow(data)==0){
    p<-barplot(ego, showCategory=20);
    print(p);
    edox<-setReadable(ego,OrgDb=org.Rn.eg.db,keyType = "ENTREZID")
    write.table(edox,paste0("AdZ4G1_vs_AdVector--up_in_",i,"_KEGG_",type,".csv"))
  }
  dev.off()
  # Gene-Concept Network
  pdf(paste0("AdZ4G1_vs_AdVector--up_in_",i,"_Gene-Concept-Network.pdf",sep=""))
  input<-AdZ4G1_vs_AdVector$log2FoldChange
  names(input)<-rownames(AdZ4G1_vs_AdVector)
  p1 <- cnetplot(edox, foldChange=input, cex_label_category=0.8)
  p3 <- cnetplot(edox, foldChange=input, circular = TRUE, colorEdge = TRUE) 
  print(p1);
  print(p3);
  dev.off()
}
#GSEA 
library(enrichplot)
gene=bitr(rownames(AdZ4G1_vs_AdVector),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Rn.eg.db") 

gene_df <- data.frame(logFC=AdZ4G1_vs_AdVector$log2FoldChange, 
                      SYMBOL = rownames(AdZ4G1_vs_AdVector)) 
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL
geneList=sort(geneList,decreasing = T) 
kegmt<-read.gmt("D:/project/HMMR/HMMR-OE-P1/rattus_norvegicus_brown-rat_gmt2.gmt") 
names(geneList)<-toupper(names(geneList))
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.45,
           TERM2GENE = kegmt) 
write.csv(KEGG,"AdZ4G1_vs_AdVector_GSVA.csv")
pdf("AdZ4G1_vs_AdVector_GSVA-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="p.adjust",split=".sign")+facet_wrap(~.sign,scales = "free") 
dev.off()
pdf("AdZ4G1_vs_AdVector_GSVA.pdf",width = 10,height = 10)
for (i in 1:200){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) 
dev.off()
library(STRINGdb) 
library(tidyverse) 
library(clusterProfiler) 
library(org.Rn.eg.db) 
library(igraph) 
library(ggraph)
string_db <- STRINGdb$new( version="11.5",species=10116)   
STRINGdb$methods() 
for (i in c("AdZ4G1","AdVector")){
  print(i);
  gene_in_cluster<- rownames(res_data[res_data$change==i,])
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"),OrgDb = org.Rn.eg.db)
  #GO
  pdf(paste0("AdZ4G1_vs_AdVector--up_in_",i,"-Stringdb-network.pdf",sep=""))
  {data_mapped <- gene.df %>% string_db$map(my_data_frame_id_col_names = "ENTREZID",removeUnmappedRows = TRUE)  
    string_db$plot_network( data_mapped$STRING_id )
    data_links <- data_mapped$STRING_id %>% string_db$get_interactions() 
    links <- data_links %>%
      mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
      mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
      dplyr::select(from, to , last_col()) %>% 
      dplyr::rename(weight = combined_score)
    nodes <- links %>% {data.frame(gene = c(.$from, .$to)) } %>% distinct()
    net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
    igraph::V(net)$deg <- igraph::degree(net) # 每个节点连接的节点数
    igraph::V(net)$size <- igraph::degree(net)/5 #
    igraph::E(net)$width <- igraph::E(net)$weight/10
    ggraph(net,layout = "kk")+
      geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
      geom_node_point(aes(size=size), color="orange", alpha=0.7)+
      geom_node_text(aes(filter=deg>1, label=name), size = 5, repel = T)+
      scale_edge_width(range = c(0.2,1))+
      scale_size_continuous(range = c(1,10) )+
      guides(size=F)+
      theme_graph()
    dev.off()
  }
}


#3D PCA#
library("FactoMineR")
library("factoextra")
library("scatterplot3d")
library("gmodels")
pca_count<- limma::removeBatchEffect(assay(vsd), c(colData$batch))
#pca_count<- normalized_counts
pca.info <- fast.prcomp(pca_count)
head(pca.info$rotation)
head(pca.info)
pca.data <- data.frame(sample =rownames(pca.info$rotation),
                       condition=condition,
                       pca.info$rotation)

library(ggsci)
library("scales")
colors=pal_npg("nrc")(10)
show_col(pal_npg("nrc")(10))
colors_pal<-colors[c(3,2,5,1)]
colors <- colors_pal[as.factor(pca.data$condition)]

pVar <- pca.info$sdev^2/sum(pca.info$sdev^2)
pVar = round(pVar,digits = 3)

paste0("PC1 (",as.character(pVar[1] * 100 ),"%)")
paste0("PC2 (",as.character(pVar[2] * 100 ),"%)")
paste0("PC3 (",as.character(pVar[3] * 100 ),"%)")

str(pca.data)
s3d <- scatterplot3d(pca.data[,c("PC1","PC2","PC3")],pch=15,color = colors,mar=c(5,5,5,5),
                     angle = 60, type="p",cex.symbols = 1,
                     main = "3D PCA plot",
                     xlab="PC1 ",
                     ylab = "PC2",
                     zlab = "PC3 ") 
legend("topright", legend = c("AdG1","AdVector","AdZ4","AdZ4G1"),
       col =colors_pal, pch = 15, bty="n",
       inset = -0.2,xpd = TRUE, horiz = FALSE)
