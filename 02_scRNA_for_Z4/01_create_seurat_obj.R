## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
library(DoubletFinder)
##Load datasets and create Seurat objects with the raw (non-normalized data).
AR_G1   <-  CreateSeuratObject(counts = Read10X(data.dir = "/md01/nieyg/project/Maternal_scRNA/00_data/1_cellranger/AR-G1-10SN/outs/filtered_feature_bc_matrix"), 
  project="AR_G1", assay = "RNA")
AR_GFP <-  CreateSeuratObject(counts = Read10X(data.dir = "/md01/nieyg/project/Maternal_scRNA/00_data/1_cellranger/AR-GFP-10SN/outs/filtered_feature_bc_matrix"), 
  project="AR_GFP", assay = "RNA")
AR_Z4 <-  CreateSeuratObject(counts = Read10X(data.dir = "/md01/nieyg/project/Maternal_scRNA/00_data/1_cellranger/AR-Z4-10SN/outs/filtered_feature_bc_matrix"), 
  project="AR_Z4", assay = "RNA")
Sham <-  CreateSeuratObject(counts = Read10X(data.dir = "/md01/nieyg/project/Maternal_scRNA/00_data/1_cellranger/Sham-10SN/outs/filtered_feature_bc_matrix"), 
  project="Sham", assay = "RNA")
#####添加sample信息####

##############################
#######pre-processing#########
#######      1_QC      #########
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#注意参考基因组里面线粒体相关基因的名称
AR_G1[["percent.mt"]] <- PercentageFeatureSet(AR_G1, pattern = "^mt-")
AR_GFP[["percent.mt"]] <- PercentageFeatureSet(AR_GFP, pattern = "^mt-")
AR_Z4[["percent.mt"]] <- PercentageFeatureSet(AR_Z4, pattern = "^mt-")
Sham[["percent.mt"]] <- PercentageFeatureSet(Sham, pattern = "^mt-")

#devtools::install_github("thomasp85/patchwork") 

# Visualize QC metrics as a violin plot
pdf("./01_QC/sample_qc_plot_nopoint.pdf")
VlnPlot(AR_GFP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
VlnPlot(AR_Z4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
VlnPlot(AR_G1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
VlnPlot(Sham, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("./01_QC/samplel_qc_plot.pdf")
plot3 <- FeatureScatter(Sham, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(Sham, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
plot3 <- FeatureScatter(AR_GFP, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(AR_GFP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
plot1 <- FeatureScatter(AR_G1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AR_G1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot3 <- FeatureScatter(AR_Z4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(AR_Z4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
dev.off()
######################################################################################
# 2_DOUBLETFINDER
######################################################################################

# AR_Z4MALIZE for DoubletFinder
library(gridExtra)
library(DoubletFinder)
sample <- c("AR_GFP","AR_G1","AR_Z4","Sham")
objList <- list(AR_GFP,AR_G1,AR_Z4,Sham)

for (i in seq_len(length(objList))) {
  # Normalizing
  dataset.name <- sample[i]
  sratDecontx    <- SCTransform(objList[[i]], verbose = F)
  # Run cluster analysis
  sratDecontx    <- RunPCA(sratDecontx, npcs = 50, verbose = F)
  sratDecontx    <- RunUMAP(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- RunTSNE(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- FindNeighbors(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- FindClusters(sratDecontx, verbose = T)
  # PLOT: umap/tsne
  p1 <- DimPlot(sratDecontx, reduction = "umap", label = TRUE)
  p2 <- DimPlot(sratDecontx, reduction = "tsne", label = TRUE)
  g = arrangeGrob(p1,p2, ncol = 2)
  filename <- paste0("./01_QC/1b_filtering_UMAP-tSNE_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,g)
  rm(g)
  objList[[i]] <- sratDecontx
    }


objList2<-list()
for (i in seq_len(length(objList))) {
  dataset.name<-sample[i]
  print(dataset.name)
  sratDecontx <- objList[[i]]
  # Compute expected doublet rate
  cellN=nrow(sratDecontx@meta.data)
  expDoubletRate = (cellN*0.0008 + 0.0527)*0.01
  normalizationMethod='SCTransform'
  sweep.res.list_scData <- paramSweep_v3(sratDecontx, 
                                         PCs = 1:50, 
                                         sct = TRUE,num.cores = 4) #num.cores = 4
  sweep.stats_scData <- summarizeSweep(sweep.res.list_scData, GT = FALSE)
  bcmvn_scData <- find.pK(sweep.stats_scData)
  bcmvn_scData$pK <- as.numeric(as.character(bcmvn_scData$pK))
  pK1=bcmvn_scData$pK[bcmvn_scData$BCmetric==max(bcmvn_scData$BCmetric)]
  print(head(pK1))
  # PLOT: pK selection
  p1=ggplot(data=bcmvn_scData, 
            aes(x=pK, y=BCmetric, group=2)) +
    geom_line(color="blue")+
    geom_point()+
    geom_vline(xintercept=pK1, linetype="dashed", color = "red")+
    labs(title="pK Selection",x="pK", y = "BCmvn")+
    theme_classic()
  filename <- paste0("./01_QC/2a_doubletfinder_pkselection_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p1)
  # More doublet finder
  pK1=as.numeric(as.character( pK1 ))
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- sratDecontx@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)   
  nExp_poi <- round(expDoubletRate*nrow(sratDecontx@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  sratDecontx <- doubletFinder_v3( sratDecontx, PCs = sratDecontx@commands$RunUMAP.SCT.pca$dims,
                                   pN = 0.25, pK = pK1, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  sratDecontx@meta.data$DoubletFinder =  sratDecontx@meta.data[,grep('DF.classifications', colnames( sratDecontx@meta.data))]
  # PLOT: Doublet Finder graphs
  p2 <- FeatureScatter(sratDecontx, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'DoubletFinder')
  p3 <- FeatureScatter(sratDecontx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'DoubletFinder')
  g = arrangeGrob(p2,p3, ncol = 2)
  filename <- paste0("./01_QC/2b_doubletfinder_ScatterPlots_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,g)
  rm(filename, g)
  # PLOT: Violin Plots
  p4 <- VlnPlot(sratDecontx, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), group.by = 'DoubletFinder', pt.size = 0)
  filename <- paste0("./01_QC/2c_doubletfinder_DoubletFinder-ViolinPlots_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p4)
  objList2[[i]] <- sratDecontx
}

######################################################################################
# 3_FEATURE FILTERING
######################################################################################

objList<- objList2
# Subset identified singlets
for (i in seq_len(length(objList))) {
    objList[[i]] <- subset(objList[[i]], subset = DoubletFinder == "Singlet")
}


library(ggplot2)

for (i in seq_len(length(objList))) {
    objList[[i]] <- subset(objList[[i]],nFeature_RNA >200 & nFeature_RNA < 8000 & percent.mt < 20 )
}
#objList <- list(AR_GFP,AR_G1,AR_Z4)
AR_GFP<- subset(AR_GFP,cells=colnames(objList[[1]]))
AR_G1<- subset(AR_G1,cells=colnames(objList[[2]]))
AR_Z4<- subset(AR_Z4,cells=colnames(objList[[3]]))
Sham<- subset(Sham,cells=colnames(objList[[4]]))

# Simply merge Seurat objects
merged_obj <- merge(x=AR_Z4,y=c(AR_GFP,Sham),add.cell.ids = c("AR_Z4","AR_GFP","Sham"),project = "Maternal")
Idents(merged_obj) <- merged_obj$orig.ident
#split the combined object into a list, with each dataset as an element
Maternal.list <- SplitObject(merged_obj,split.by = "ident")

pdf("./01_QC/AR_Z4_normalization_plot.pdf")
Maternal.list[[1]] <- FindVariableFeatures(Maternal.list[[1]], verbose = FALSE)
top10 <- head(VariableFeatures(Maternal.list[[1]]),10)
plot1 <- VariableFeaturePlot(Maternal.list[[1]])
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
dev.off();

Maternal.anchors <- FindIntegrationAnchors(object.list = Maternal.list,anchor.features = 2000,dims = 1:30)
Maternal.integrated <- IntegrateData(anchorset = Maternal.anchors, dims = 1:30,features.to.integrate = rownames(Maternal.list[[1]]))
head(Maternal.integrated@meta.data)
DefaultAssay(Maternal.integrated) <- "integrated"

# scale and center features in the dataset
Maternal.integrated <- ScaleData(Maternal.integrated, features =rownames(Maternal.integrated),verbose = FALSE)
# Perform linear dimensional reduction
Maternal.integrated <- RunPCA(Maternal.integrated, npcs = 50, verbose = FALSE)
# Determine the ‘dimensionality’ of the dataset
# JackStraw 
#Maternal.integrated <- JackStraw(Maternal.integrated, num.replicate = 100, dims =50)
Maternal.integrated  <- ScoreJackStraw(Maternal.integrated, dims = 1:50)
pdf("./03_Z4/01_All_celltype/seclect_pca.pdf")
#JackStrawPlot(Maternal.integrated, dims = 1:50)
# ‘Elbow plot’
ElbowPlot(Maternal.integrated,ndims=50)
dev.off()
library(clustree)
resolutions <- seq(0.2, 1, 0.1)

DefaultAssay(Maternal.integrated) <- "integrated"
Maternal.integrated <- FindNeighbors(object = Maternal.integrated, dims = 1:25)
Maternal.integrated <- FindClusters(object = Maternal.integrated, resolution = resolutions)
Maternal.integrated <- RunUMAP(object = Maternal.integrated, dims = 1:25)

pdf("./03_Z4/01_All_celltype/resolutions_select_Z4_all_celltype.pdf",width=15,height=15)
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
Maternal.integrated$seurat_clusters<- Maternal.integrated$integrated_snn_res.0.3
Idents(Maternal.integrated)<- Maternal.integrated$seurat_clusters

pdf("./03_Z4/01_All_celltype/Unsupervised_cluster_Umap.pdf")
DimPlot(object = Maternal.integrated, group.by = "orig.ident",label = T)
DimPlot(object = Maternal.integrated, group.by = "seurat_clusters",label = T)
dev.off()
DefaultAssay(Maternal.integrated) <- "RNA"
saveRDS(Maternal.integrated, file = "./03_Z4/01_All_celltype/Maternal_integrated_Z4.rds")




