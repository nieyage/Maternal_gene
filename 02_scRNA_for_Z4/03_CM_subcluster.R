## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
set.seed(1234)
combined<- readRDS("./01_All_celltype/Maternal_Z4_integrated_annotated.rds")

CM<- subset(combined,ident="CM")
Idents(CM)<- CM$orig.ident
Z4.list <- SplitObject(CM,split.by = "ident")
for (i in 1:3){
  Z4.list[[i]] <- FindVariableFeatures(Z4.list[[i]], verbose = FALSE)
}

Z4.anchors <- FindIntegrationAnchors(object.list = Z4.list,anchor.features = 2000,dims = 1:30)
Z4.integrated <- IntegrateData(anchorset = Z4.anchors, dims = 1:30,new.assay.name = "integrated_CM")
head(Z4.integrated@meta.data)
DefaultAssay(Z4.integrated) <- "integrated_CM"

# scale and center features in the dataset
Z4.integrated <- ScaleData(Z4.integrated, features =rownames(Z4.integrated),verbose = FALSE)
# Perform linear dimensional reduction
Z4.integrated <- RunPCA(Z4.integrated, npcs = 50, verbose = FALSE)
# Determine the ‘dimensionality’ of the dataset

#Z4.integrated <- JackStraw(Z4.integrated, num.replicate = 100, dims =50)
Z4.integrated  <- ScoreJackStraw(Z4.integrated, dims = 1:50)
pdf("./02_CM/seclect_pca.pdf")
#JackStrawPlot(Z4.integrated, dims = 1:50)
# ‘Elbow plot’
ElbowPlot(Z4.integrated,ndims=50)
dev.off()

DefaultAssay(Z4.integrated) <- "integrated_CM"
resolutions <- seq(0.2, 2, 0.2)
CM<- Z4.integrated
CM <- FindNeighbors(object = CM, dims = 1:20)
CM <- FindClusters(object = CM, resolution = resolutions)
CM <- RunUMAP(object = CM, dims = 1:20)

pdf("./02_CM/resolutions_select_G1_CM_subcluster.pdf",width=15,height=15)
clustree(CM, prefix = "integrated_CM_snn_res.")
plotlist <- lapply(resolutions, function(x){
    cols <- myUmapcolors
    p <- DimPlot(CM, group.by = glue::glue("integrated_CM_snn_res.{x}"), label = TRUE,
             reduction = "umap", shuffle = TRUE) +
    scale_color_manual(values = cols) +
    xlab("UMAP1") + ylab("UMAP2")
    p
})
p <- patchwork::wrap_plots(plotlist, ncol = 3)
p
dev.off()
CM$seurat_clusters<- CM$integrated_CM_snn_res.0.2
Idents(CM)<- CM$seurat_clusters

pdf("./02_CM/CM_Unsupervised_cluster_Umap.pdf")
DimPlot(object = CM, group.by = "orig.ident")
DimPlot(object = CM,cols=myUmapcolors, group.by = "seurat_clusters")
dev.off()

# make the tree 
# make the trans dist tree 
object <- CM
embeddings <- Embeddings(object = object, reduction = "pca")[,1:20]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
library(ggtree);

pdf("./02_CM/CM_cluster-tree-cosine.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()
DefaultAssay(CM)<-"RNA"
pdf('./02_CM/CM_markers_VlnPlot.pdf', width=12, height=8)
VlnPlot(CM,cols=myUmapcolors, features = c("Tnnt2","Myl2","Actc1","Myh6"), ncol = 1, pt.size = 0)
VlnPlot(CM,cols=myUmapcolors, y.max=20,features = c("Tnnt2","Myl2","Actc1","Myh6"), ncol = 1, pt.size = 0)
VlnPlot(CM,cols=myUmapcolors, features = c("nCount_RNA","nFeature_RNA","percent.mt"), ncol = 1, pt.size = 0)
dev.off()

# plot the features 
pdf("./02_CM/CM_subcluster_umap_features.pdf",width=10,height=10)
FeaturePlot(
  object = CM,
  features = c("nCount_RNA","nFeature_RNA","percent.mt"),
  order = FALSE,
  reduction='umap' ,
  cols = c("grey", "darkred"),
  ncol = 2
) & NoLegend()
dev.off()

pdf("./02_CM/CM_subcluster_umap_features_Zscan4c.pdf",width=10,height=10)
FeaturePlot(
  object = CM,
  features = c("Zscan4c","Tnnt2","Myl2","Actc1"),
  order = TRUE,
  reduction='umap' ,
  cols = c("grey", "darkred"),
  ncol = 2
) & NoLegend()
dev.off()

combined<- CM
Idents(combined)<-combined$seurat_clusters
#####further annotation########
combined <- RenameIdents(
  object = combined,
  '0' = 'CM1',
  '1' = 'CM2',
  '2' = 'CM3',
  '3' = 'CM4',
  '4' = 'CM3',
  '5' = 'CM5',
  '6' = 'CM6',
  '7' = 'CM7',
  '8' = 'CM8'
  )
combined@meta.data$subtype<-Idents(combined)
table(combined$subtype,combined$orig.ident)
Idents(combined)<-combined$subtype;
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

pdf("./02_CM/CM_subtype_UMAP.pdf",width=5,height=5)
DimPlot(combined,  repel = TRUE, cols=myUmapcolors[1:8], reduction = "umap",group.by = "subtype")
dev.off()


cell.prop<-as.data.frame(prop.table(table(Idents(combined), combined$orig.ident)))
colnames(cell.prop)<-c("subtype","sample","proportion")
cell.prop$sample<- factor(cell.prop$sample,levels=c("Sham","AR_GFP","AR_Z4"))

pdf("./02_CM/CM_subtype_proportion.pdf",width=4,height=5)
ggplot(cell.prop,aes(sample,proportion,fill=subtype))+
geom_bar(stat="identity",position="fill")+
ggtitle("")+
theme_bw()+
theme(axis.ticks.length=unit(0.5,'cm'))+
guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values = myUmapcolors[1:8]) 
dev.off()

# Save rds have annotation information 
DefaultAssay(combined) <- "RNA"
saveRDS(combined,"./02_CM/CM_subtype_annotated.rds")

combined<- readRDS("./02_CM/CM_subtype_annotated.rds")
combined<- SCTransform(combined)
DefaultAssay(combined)<- "SCT"
combined<- ScaleData(combined)

markers <- FindAllMarkers(combined, only.pos = TRUE,logfc.threshold = 0.25)
#markers<- markers[markers$p_val_adj<0.05&markers$avg_log2FC>1,]

write.csv(markers,"./02_CM/CM_DEG_markers-log2FC1adjP005.csv")

top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC);

blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")
pdf("./02_CM/CM_markers-log2FC1adjP005_top_markers-DEG_heatmap.pdf",width=10,height=10)
DoHeatmap(object = combined,features=markers$gene,label=T,size = 2,group.by = "subtype") 
DoHeatmap(object = combined,features=top10$gene,label=T, group.colors =c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494","#B3B3B3" ),
  disp.min = -2,disp.max = 2,size = 2,group.by = "subtype") 
DoHeatmap(object = combined,features=markers$gene,label=T, group.colors =c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494","#B3B3B3" ),
  disp.min = -2,disp.max = 2,size = 2,group.by = "subtype") 
dev.off();

Idents(combined)<- combined$subtype
CM_sub<- subset(combined,ident=c("CM1","CM2","CM3"))
DefaultAssay(CM_sub)<- "integrated_CM"
CM_sub <- FindNeighbors(object = CM_sub, dims = 1:20)
CM_sub <- FindClusters(object = CM_sub, resolution = 1)
CM_sub <- RunUMAP(object = CM_sub, dims = 1:20)

DefaultAssay(CM_sub)<- "RNA"
pdf("./02_CM/CM_subcluster2_umap_features_Zscan4c.pdf",width=10,height=10)
FeaturePlot(
  object = CM_sub,
  features = c("Zscan4c","Tnnt2","Myl2","Actc1"),
  order = TRUE,
  reduction='umap' ,
  cols = c("grey", "darkred"),
  ncol = 2
) & NoLegend()
dev.off()

























# GO and KEGG for CM subtype
CM1<-markers[markers$cluster=="CM1",7]
CM2<-markers[markers$cluster=="CM2",7]
CM3<-markers[markers$cluster=="CM3",7]
CM4<-markers[markers$cluster=="CM4",7]
CM5<-markers[markers$cluster=="CM5",7]

library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

CM_list<- list(CM1,CM2,CM3,CM4,CM5)
subtype<- c("CM1","CM2","CM3","CM4","CM5")
all_ego<- data.frame()
pdf("./02_CM/CM_subtype_DEG_GO_BP.pdf")
for(i in 1:length(subtype)){
	gene<- CM_list[[i]];
	gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.5,
                readable = TRUE)
    p<- barplot(ego, showCategory=20)
    print(p)
    write.csv(ego,paste0(subtype[i],"-GO-BP.csv",sep=""))
    ego<- as.data.frame(ego)
    ego$celltype<- subtype[i];
    all_ego<- rbind(all_ego,ego);

}
dev.off()

write.csv(all_ego,"./02_CM/Allsubtype_GO.csv")

top_term <- read.csv("/data/R02/nieyg/project/Z4_scRNA/01_all_celltype/02_CM/Allsubtype_GO-selected.csv") 
library(ggplot2)
top_term$Description<- factor(top_term$Description,levels=rev(top_term$Description))
pdf("/data/R02/nieyg/project/Z4_scRNA/01_all_celltype/02_CM/Allcelltype_GO.pdf",width=14,height=12)
p <- ggplot(top_term,aes(y=Count,x=Description,fill=pvalue)) + 
      geom_bar(stat="identity",position = "dodge") +
      facet_grid(celltype~.,scales = "free",space = "free") + 
      coord_flip() + 
      theme_bw() +
      scale_fill_gradient(low = '#FF0000', high = '#1202FF')+
      theme(plot.title = element_text(hjust = 0.5),
            strip.text.y = element_text(size = 14),
            legend.position="right",
            legend.title = element_text(size=18),
            legend.text = element_text(size=14),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=18),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
p
dev.off()



library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

CM_list<- list(CM1,CM2,CM3,CM4,CM5)
subtype<- c("CM1","CM2","CM3","CM4","CM5")
all_ego<- data.frame()
pdf("./CM_subtype_DEG_KEGG.pdf")
for(i in 1:length(subtype)){
    gene<- CM_list[[i]];
  gene.df <- bitr(gene, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  ego <- enrichKEGG(gene.df$ENTREZID,
                    keyType = "kegg",
                    organism  = 'hsa',
                    pvalueCutoff  = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 0.5)
  p<- barplot(ego, showCategory=20)
  print(p)
  ego<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
  write.csv(ego,paste0("./CMsubtype-",subtype[i],"-KEGG.csv",sep=""))
    ego<- as.data.frame(ego)
    ego$celltype<- subtype[i];
    all_ego<- rbind(all_ego,ego);

}
dev.off()
write.csv(all_ego,"./CMsubtype_KEGG.csv")




