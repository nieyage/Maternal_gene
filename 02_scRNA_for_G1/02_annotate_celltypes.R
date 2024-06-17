## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
set.seed(1234)
Maternal.integrated<- readRDS("./01_All_celltype/Maternal_integrated_G1.rds")
DefaultAssay(Maternal.integrated) <- "integrated"
library(clustree)
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

resolutions <- seq(0.2, 1, 0.1)
Maternal.integrated <- FindNeighbors(object = Maternal.integrated, dims = 1:25)
Maternal.integrated <- FindClusters(object = Maternal.integrated, resolution = resolutions)
Maternal.integrated <- RunUMAP(object = Maternal.integrated, dims = 1:25)
table(Maternal.integrated@active.ident)

pdf("./01_All_celltype/resolutions_select_G1_all_celltype.pdf",width=15,height=15)
clustree(Maternal.integrated, prefix = "integrated_snn_res.")
plotlist <- lapply(resolutions, function(x){
    cols <- myUmapcolors
    p <- DimPlot(Maternal.integrated, group.by = glue::glue("integrated_snn_res.{x}"), label = TRUE,
             reduction = "umap", shuffle = TRUE) +
    scale_color_manual(values = cols) +
    xlab("UMAP1") + ylab("UMAP2")
    p
})
p <- patchwork::wrap_plots(plotlist, ncol = 3)
p
dev.off()
Maternal.integrated$seurat_clusters<- Maternal.integrated$integrated_snn_res.0.4
Idents(Maternal.integrated)<- Maternal.integrated$seurat_clusters

pdf("./01_All_celltype/Unsupervised_cluster_Umap.pdf")
DimPlot(object = Maternal.integrated, group.by = "orig.ident",label = T)
DimPlot(object = Maternal.integrated, group.by = "seurat_clusters",label = T)
dev.off()

table(Maternal.integrated@active.ident)
#cluster annotation
Idents(Maternal.integrated)<- Maternal.integrated$seurat_clusters
DefaultAssay(Maternal.integrated) <- "RNA"
markerGenes  <- c(
"Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1", # cardiomyocyte
"Ccl3","Clec4a1","Rgs1","Fcgr1","Cd14","Csf1r","Cd163","Cd68","Itgam","Lgals3","Mrc1", # Macrophages
"Cd3e","Cd3d","Cd8a","Cd8b1","Igfbp4","Lat","Itk","Cd3g", # T cell
"Cd22","Cd79a","Cd79b","Mzb1","Ly6d","Ms4a1", # B cell
"Cd74","Cd83","Cd86","Flt3","Cd209a","Ccl17", # DC
"Plac8", # monocyte
"S100a9", "S100a8",#Granulocy
"Nkg7","Gzma", #NK
"Chl1","Kcna6", # Glial
"Myct1","Cdh5","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1","Tie1","Fabp4","Esam","Kdr","Tek", # endothelial
"Dpt","Col3a1","Mmp2","Col1a2","Fstl1","Gsn","Thy1","Pdgfra", # fibroblast
"Abcc9","Rgs4","Ano1","Acta2","Myh11", # smoothmuscle
"Msln","Krt19","Plet1","Prr15","Slc9a3r1","Lrrn4","Slc39a8","Krt8","Anxa8",# epicardial
"Gfpt2","Hk2","Ddr2","Slc1a7","Adamtsl3","Layn","Cd248","Mcam","Fgfr1" # Pericyte
)

label<- c(rep("CM",7),rep("MP",11),rep("T",8),rep("B",6),
    rep("DC",6),rep("monocyte",1),rep("Granulocy",2),
    rep("NK",2),rep("Glial",2),rep("EC",13),
  rep("FB",8),rep("SMC",5),rep("Epi",9),rep("Pericyte",9))
label<- label[which(markerGenes%in% rownames(Maternal.integrated))]
markerGenes<- markerGenes[which(markerGenes%in% rownames(Maternal.integrated))]

DefaultAssay(Maternal.integrated)<-"RNA"
pdf("./01_All_celltype/G1_cluster-annotation-all_celltype.pdf",width=20,height=8)
p<-DotPlot(Maternal.integrated, features = markerGenes,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
p&scale_x_discrete(labels=label)
dev.off()


DefaultAssay(Maternal.integrated) <- "RNA"
markerGenes  <- c(
"Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1", # cardiomyocyte
"Ccl3","Clec4a1","Rgs1","Fcgr1","Cd14","Csf1r","Cd163","Cd68","Itgam","Lgals3","Mrc1", # Macrophages
"Cd3e","Cd3d","Cd8a","Cd8b1","Igfbp4","Lat","Itk","Cd3g", # T cell
"Cd22","Cd79a","Cd79b","Mzb1","Ly6d","Ms4a1", # B cell
"Cd74","Cd83","Cd86","Flt3","Cd209a","Ccl17", # DC
"Plac8", # monocyte
"S100a9", "S100a8","Wfdc21","Ly6c2",#Granulocyte
"Nkg7","Gzma", #NK
"Chl1","Kcna6", # Glial
"Myct1","Cdh5","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1","Tie1","Fabp4","Esam","Kdr","Tek", # endothelial
"Dpt","Col3a1","Mmp2","Col1a2","Fstl1","Gsn","Thy1","Pdgfra", # fibroblast
"Abcc9","Rgs4","Ano1","Acta2","Myh11", # smoothmuscle
"Msln","Krt19","Plet1","Prr15","Slc9a3r1","Lrrn4","Slc39a8","Krt8","Anxa8",# epicardial
"Gfpt2","Hk2","Ddr2","Slc1a7","Adamtsl3","Layn","Cd248","Mcam","Fgfr1","Rgs5","Kcnj8","Vtn", # Pericyte
"C1qa",# myeloid
"Il7r","Cd40lg",# Lymphoid
"Plp1","Nrxn1","Nrxn3",# Neuronal
"Gpam","Fasn","Lep",#Adipocytes
"Wt1","Bnc1"# Mesothelial
)

label<- c(rep("CM",8),rep("MP",11),rep("T",8),rep("B",6),
    rep("DC",6),rep("monocyte",1),rep("Granulocy",4),
    rep("NK",2),rep("Glial",2),rep("EC",13),
  rep("FB",8),rep("SMC",5),rep("Epi",9),rep("Pericyte",12),
  rep("myeloid",1),rep("Lymphoid",2),rep("Neuronal",3),rep("Adipocytes",3),rep("Mesothelial",2))
label<- label[which(markerGenes%in% rownames(Maternal.integrated))]
markerGenes<- markerGenes[which(markerGenes%in% rownames(Maternal.integrated))]

DefaultAssay(Maternal.integrated)<-"RNA"
pdf("./01_All_celltype/G1_cluster-annotation-all_celltype2.pdf",width=16,height=8)
p<-DotPlot(Maternal.integrated, features = markerGenes,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
p&scale_x_discrete(labels=label)
dev.off()

# make the tree 

# make the trans dist tree 
object <- Maternal.integrated
embeddings <- Embeddings(object = object, reduction = "pca")[,1:50]
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

pdf("./01_All_celltype/G1_Maternal_cluster-tree-cosine.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()

combined<- Maternal.integrated
Idents(combined)<-combined$seurat_clusters
#####further annotation########
combined <- RenameIdents(
  object = combined,
  '0' = 'FB',
  '1' = 'CM',
  '2' = 'EC',
  '3' = 'MP',
  '4' = 'SMC',
  '5' = 'FB',
  '6' = 'FB',
  '7' = 'EC',
  '8' = 'EC',
  '9' = 'EC',
  '10' = 'Epi',
  '11' = 'SMC',
  '12' = 'B',
  '13' = 'EC',
  '14' = 'MP',
  '15' = 'MP',
  '16' = 'unknown',
  '17' = 'T',
  '18' = 'SMC',
  '19' = 'Glial',
  '20' = 'MP',
  '21' = 'unknown'
  )
combined@meta.data$Annotation<-Idents(combined)
table(combined$Annotation,combined$orig.ident)
Idents(combined)<- combined$Annotation
combined_sub<- subset(combined,ident=c("CM",'EC','FB','SMC','Epi','MP','T','B','Glial'))

combined<- combined_sub
combined$Annotation<-factor(combined$Annotation,levels=c("CM",'EC','FB','SMC','Epi','MP','T','B','Glial'))
Idents(combined)<-combined$Annotation;
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

pdf("./01_All_celltype/G1_Annotated_allcelltype_UMAP.pdf",width=18,height=5)
DimPlot(combined, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap",group.by = c("orig.ident","seurat_clusters","Annotation"))
DimPlot(combined, label = T, repel = TRUE, reduction = "umap",group.by = c("orig.ident","seurat_clusters","Annotation"))
dev.off()

cell.prop<-as.data.frame(prop.table(table(Idents(combined), combined$orig.ident)))
colnames(cell.prop)<-c("celltype","sample","proportion")
cell.prop$sample<- factor(cell.prop$sample,levels=c("Sham","AR_GFP","AR_G1"))
cell.prop$celltype<- factor(cell.prop$celltype,levels=c("CM",'EC','FB','SMC','Epi','MP','T','B','Glial'))

pdf("./01_All_celltype/G1_Annotated_allcelltype_proportion.pdf",width=4,height=5)
ggplot(cell.prop,aes(sample,proportion,fill=celltype))+
geom_bar(stat="identity",position="fill")+
ggtitle("")+
theme_bw()+
theme(axis.ticks.length=unit(0.5,'cm'))+
guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values = myUmapcolors) 
dev.off()

# Save rds have annotation information 
DefaultAssay(combined) <- "RNA"
saveRDS(combined,"./01_All_celltype/Maternal_G1_integrated_annotated.rds")

combined<- readRDS("./01_All_celltype/Maternal_G1_integrated_annotated.rds")

# plot the marker gene among different celltypes 
# dotplot 
 
DefaultAssay(combined) <- "RNA"
Idents(combined)<- combined$Annotation
markerGenes  <- c(
"Tnnt2","Myl2","Actc1","Myh6","Atp2a2", # cardiomyocyte
"Myct1","Cdh5","Flt1","Tie1","Fabp4",# endothelial
"Col3a1","Mmp2","Dpt","Thy1","Gsn", # fibroblast
"Abcc9","Rgs4","Ano1","Acta2","Myh11", # smoothmuscle
"Msln","Krt19","Slc9a3r1","Lrrn4","Slc39a8",# epicardial
"Csf1r","Cd163","Itgam","Lgals3","Mrc1", # Macrophages
"Cd3e","Cd3d","Cd8a","Itk","Cd3g", # T cell
"Cd22","Cd79a","Cd79b","Ly6d","Ms4a1", # B cell
"Fgfr1", # Pericyte
"Chl1","Kcna6"# Glial
)

label<- c(rep("CM",5),rep("EC",5),
  rep("FB",5),rep("SMC",5),rep("Epi",5),
  rep("MP",5),rep("T",5),rep("B",5),
  rep("Pericyte",1),rep("Glial",2))
label<- label[which(markerGenes%in% rownames(Maternal.integrated))]
markerGenes<- markerGenes[which(markerGenes%in% rownames(Maternal.integrated))]

DefaultAssay(combined) <- "RNA"
Idents(combined)<- combined$Annotation
pdf("./01_All_celltype/G1_cluster-annotated-dotplot.pdf",width=10,height=5)
p<-DotPlot(combined, features = markerGenes,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
p&scale_x_discrete(labels=label)
dev.off()

# heatmap for specific gene top 10 
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
write.csv(markers,"./01_All_celltype/G1_FindAllMarkers_celltypes.csv")
combined<-ScaleData(combined,features=rownames(combined));
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")

pdf("./01_All_celltype/G1_Allcelltype_top_markers-DEG_heatmap.pdf",width=10,height=10)
DoHeatmap(object = combined,features=top10$gene,label=T, group.colors = myUmapcolors[1:10],
  disp.min = -2,disp.max = 2,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = combined,features=top10$gene,label=T, group.colors = myUmapcolors[1:10],
  disp.min = -2,disp.max = 2,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = solarExtra)
DoHeatmap(object = combined,features=top100$gene,label=T, group.colors = myUmapcolors[1:10],
  disp.min = -2,disp.max = 2,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = combined,features=top100$gene,label=T, group.colors = myUmapcolors[1:10],
  disp.min = -2,disp.max = 2,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = solarExtra)
dev.off();

# GO for each celltype:
# GO for each celltype:
CM<-markers[markers$cluster=="CM",7]
EC<-markers[markers$cluster=="EC",7]
FB<-markers[markers$cluster=="FB",7]
SMC<-markers[markers$cluster=="SMC",7]
Epi<-markers[markers$cluster=="Epi",7]
MP<-markers[markers$cluster=="MP",7]
T<-markers[markers$cluster=="T",7]
B<-markers[markers$cluster=="B",7]
Glial<-markers[markers$cluster=="Glial",7]

library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
DEG_list<- list(CM,EC,FB,SMC,Epi,MP,T,B,Glial)
subtype<- c("CM",'EC','FB','SMC','Epi','MP','T','B','Glial')
all_ego<- data.frame()
pdf("./01_All_celltype/G1_All_celltype_DEG_GO_BP.pdf")
for(i in 1:length(subtype)){
  gene<- DEG_list[[i]];
  gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Mm.eg.db)
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.5,
                readable = TRUE)
    p<- barplot(ego, showCategory=20,main=subtype[i])
    print(p)
    write.csv(ego,paste0(subtype[i],"-GO-BP.csv",sep=""))
    ego<- as.data.frame(ego)
    ego$celltype<- subtype[i];
    all_ego<- rbind(all_ego,ego);
}
dev.off()
write.csv(all_ego,"./01_All_celltype/G1_All_celltype_DEG_GO_BP.csv")

# top_term <- read.csv("/data/R02/nieyg/project/Z4_scRNA/01_all_celltype/02_CM/Allsubtype_GO-selected.csv") 
library(ggplot2)
top_term <- all_ego %>% group_by(celltype) %>% top_n(n = -5, wt = p.adjust);
top_term<- top_term[!duplicated(top_term$Description),]
top_term$Description<- factor(top_term$Description,levels=rev(top_term$Description))
top_term$celltype<- factor(top_term$celltype,levels=levels(combined))
pdf("./01_All_celltype/G1_All_celltype_DEG_GO_BP_top5.pdf",width=14,height=20)
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


# Marker gene for different cell types (color on UMAP ) 
pdf('./01_All_celltype/G1_All_Marker_gene_FeaturePlot_uamp.pdf', width=5.5, height=5)
for (i in 1:length(markerGenes)){
  p<-FeaturePlot(combined,order=T, reduction = 'umap',max.cutoff = 10, features = markerGenes[i], ncol = 1)
  print(p)
}
dev.off()


