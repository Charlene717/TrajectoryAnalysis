##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages #####
  #### Basic installation ####
  Package.set <- c("tidyverse","Seurat","ggplot2","ggpmisc",
                   "stringr","magrittr","dplyr")
  ## Check whether the installation of those packages is required from basic
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)
  
  #### BiocManager installation ####
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("fgsea","AnnotationHub","ensembldb",
                   "SeuratDisk","monocle",
                   "SingleR","scRNAseq","celldex","scran")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)
  
  options(stringsAsFactors = FALSE)
  
  #### GitHub installation ####
  if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
  library(monocle)
  devtools::install_github("cole-trapnell-lab/garnett")
  devtools::install_github('cole-trapnell-lab/monocle3')
  devtools::install_github("LTLA/SingleR")
  
  library(monocle3)
  library(garnett)
  # library(SingleR)
  
  
#### Original data preprocessing #####
  load("SeuratObject_PRJCA001063.RData")
  load("H5AD_PRJCA001063_PDAC_CleanUpS_20220525.RData")
  
  seurat_meta.df <- seuratObject@meta.data
  cds_meta.df <- as.data.frame(cds@colData@listData)
  seurat_meta.df <- inner_join(seurat_meta.df, cds_meta.df)
  #seurat_meta.df[is.na(seurat_meta.df)] <- 0
  row.names(seurat_meta.df) <- seurat_meta.df[,1]
  seuratObject@meta.data <- seurat_meta.df
  DimPlot(seuratObject, reduction = "umap",group.by = "Cell_type", label = TRUE, pt.size = 0.5) + NoLegend()
  DimPlot(seuratObject, reduction = "umap",group.by = "ReCluster", label = TRUE, pt.size = 0.5) + NoLegend()
  
  sum(seuratObject@meta.data[["ReCluster"]] == "Ductal cell type 1")
  sum(seuratObject@meta.data[["ReCluster"]] == "Ductal cell type 2")
  seuratObject <- seuratObject[,!seuratObject@meta.data[["ReCluster"]] == "Ductal cell type 1"]
  seuratObject <- seuratObject[,!seuratObject@meta.data[["ReCluster"]] == "Ductal cell type 2"]
  
  seuratObject@meta.data[["ReCluster2"]] <- seuratObject@meta.data[["ReCluster"]]
  seuratObject@meta.data[["ReCluster"]] <- gsub(" ", "_", seuratObject@meta.data[["ReCluster"]])
  seuratObject@meta.data[["ReCluster"]] <- gsub("DistalCD", "MDO", seuratObject@meta.data[["ReCluster"]])
  seuratObject@meta.data[["ReCluster"]] <- gsub("CoreCD", "MDC", seuratObject@meta.data[["ReCluster"]])
  seuratObject@meta.data[["ReCluster"]] <- gsub("CDOri", "MD00", seuratObject@meta.data[["ReCluster"]])
  
  ## Modify the cell type name
  
  
  seuratObjectMono_Ori <- seuratObject
  
  library("stringr")
  rm(list=setdiff(ls(), str_subset(objects(), pattern = "seuratObject")))
  save.image("SeuratObject_CDS_PRJCA001063.RData")
  
# #### Load data #####
  load("SeuratObject_CDS_PRJCA001063.RData")
  
  
##### Current path and new folder setting* #####
  ProjectName = "TrajAna"
  Sampletype = "PDAC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }


##### Function setting #####
  ## Call function
  source("FUN_Beautify_ggplot.R") 
    

#### Plot UMAP #####
    pdf(file = paste0(Save.Path,"/",ProjectName,"_TrajectoryOri.pdf"),
        width = 20,  height = 12
    )
      Idents(seuratObject) <- seuratObject@meta.data[["Cell_type"]]
      DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
        ggtitle(paste0("CellType","  PCAi:",100,"  NNei:",20,"  MD:",0.3)) + 
        theme(plot.title = element_text(hjust = 0.5,vjust = 0)) %>% print()
      
      Idents(seuratObject) <- seuratObject@meta.data[["ReCluster"]]
      DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
        ggtitle(paste0("ReCluster","  PCAi:",100,"  NNei:",20,"  MD:",0.3)) + 
        theme(plot.title = element_text(hjust = 0.5,vjust = 0)) %>% print()
      
      FeaturePlot(seuratObject, features = c("TOP2A")) %>% BeautifyggPlot(LegPos = c(1.02, 0.15)) +
                  ggtitle(paste0("TOP2A","  PCAi:",100,"  NNei:",20,"  MD:",0.3)) + 
                  theme(plot.title = element_text(hjust = 0.5,vjust = 0)) %>% print()
    dev.off()
    
    # Re-dimension reduction
    # seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000)
    seuratObject <- FindVariableFeatures(seuratObject)
    seuratObject <- RunPCA(seuratObject,npcs = 200, features = VariableFeatures(object = seuratObject))
    ElbowPlot(seuratObject, ndims = 200)
    seuratObject <- FindNeighbors(seuratObject, dims = 1:100)
    seuratObject <- FindClusters(seuratObject, resolution = 0.5)
    
    # ## Export PDF
    # pdf(file = paste0(Save.Path,"/",ProjectName,"_Trajectory_All.pdf"),
    #     width = 20,  height = 12
    # )
    # for (i in seq(80,400,80)) {
    #   for (j in seq(0.1,0.6,0.1)) {
    #     for (k in seq(80,400,80)) {
    #       set.seed(1)
    #       seuratObject <- RunUMAP(seuratObject, dims = 1:i,n.neighbors = k, min.dist= j)
    #       seuratObject@meta.data[[paste0("UMAP_PCA",i,"_NNei",k,"_MD03",j)]] <- seuratObject@reductions[["umap"]]@cell.embeddings
    #       
    #       Idents(seuratObject) <- seuratObject@meta.data[["Cell_type"]]
    #       p <-  DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
    #         ggtitle(paste0("CellType","  PCAi:",i,"  NNei:",k,"  MD:",j)) + 
    #         theme(plot.title = element_text(hjust = 0.5,vjust = 0)) 
    #       print(p)
    #       
    #       Idents(seuratObject) <- seuratObject@meta.data[["ReCluster"]]
    #       p <-  DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
    #         ggtitle(paste0("ReCluster","  PCAi:",i,"  NNei:",k,"  MD:",j)) + 
    #         theme(plot.title = element_text(hjust = 0.5,vjust = 0)) 
    #       print(p)
    #       p <-  FeaturePlot(seuratObject, features = c("TOP2A")) %>% BeautifyggPlot(LegPos = c(1.02, 0.15)) +
    #                   ggtitle(paste0("TOP2A","  PCAi:",i,"  NNei:",k,"  MD:",j)) + 
    #                   theme(plot.title = element_text(hjust = 0.5,vjust = 0)) 
    #       print(p)
    #       
    #     }
    #   }
    # }
    # dev.off()
    
    ## Export TIFF
    for (i in seq(80,400,80)) {
      for (j in seq(0.1,0.6,0.1)) {
        for (k in seq(80,400,80)) {
          set.seed(1)
          seuratObject <- RunUMAP(seuratObject, dims = 1:i,n.neighbors = k, min.dist= j)
          seuratObject@meta.data[[paste0("UMAP_PCA",i,"_NNei",k,"_MD03",j)]] <- seuratObject@reductions[["umap"]]@cell.embeddings
          # seuratObject@reductions[["umap"]]@cell.embeddings <- seuratObject@meta.data[[paste0("UMAP_PCA",i,"_NNei",k,"_MD03",j)]]
          
          Idents(seuratObject) <- seuratObject@meta.data[["Cell_type"]]
          p <-  DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
            ggtitle(paste0("CellType","  PCAi:",i,"  NNei:",k,"  MD:",j)) + 
            theme(plot.title = element_text(hjust = 0.5,vjust = 0)) 
          tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCAi",i,"_NNei",k,"_MD",j,"_CellType.tiff"),
               width = 28, height = 20, units = "cm", res = 200)
            print(p)
          graphics.off()
          
          Idents(seuratObject) <- seuratObject@meta.data[["ReCluster"]]
          p <-  DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
            ggtitle(paste0("ReCluster","  PCAi:",i,"  NNei:",k,"  MD:",j)) + 
            theme(plot.title = element_text(hjust = 0.5,vjust = 0)) 
          tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCAi",i,"_NNei",k,"_MD",j,"_ReCluster.tiff"),
               width = 35, height = 20, units = "cm", res = 200)
            print(p)
          graphics.off()
          
          p <-  FeaturePlot(seuratObject, features = c("TOP2A")) %>% BeautifyggPlot(LegPos = c(1.02, 0.15)) +
            ggtitle(paste0("TOP2A","  PCAi:",i,"  NNei:",k,"  MD:",j)) + 
            theme(plot.title = element_text(hjust = 0.5,vjust = 0)) 
          tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCAi",i,"_NNei",k,"_MD",j,"_TOP2A.tiff"),
               width = 28, height = 20, units = "cm", res = 200)
            print(p)
          graphics.off()
          
        }
      }
    }

    

    # seuratObject <- RunUMAP(seuratObject, dims = 1:100,n.neighbors = 20, min.dist=0.05)
    # # seuratObject <- RunUMAP(seuratObject, dims = 1:100,n.neighbors = 1000, min.dist=0.1)
    # # seuratObject@meta.data[["UMAP_NNei1000"]] <- seuratObject@reductions[["umap"]]@cell.embeddings
    # # seuratObject@meta.data <- seuratObject@meta.data[,!colnames(seuratObject@meta.data)=="UMAP_NNei1000"]
    # seuratObject@meta.data[["UMAP_NNei20_MD03"]] <- seuratObject@reductions[["umap"]]@cell.embeddings
    # 
    # DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

    Idents(seuratObject) <- seuratObject@meta.data[["Cell_type"]]
    DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

    Idents(seuratObject) <- seuratObject@meta.data[["ReCluster"]]
    DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
    FeaturePlot(seuratObject, features = c("TOP2A"))

    
#### Save RData #####    
  save.image(Save.Path,"/",Version,"SeuratObject_CDS_PRJCA001063_TrajAna.RData")
    
# #### Cell type markers #####
#     ## Create cell type markers dataframe
#     # DefaultAssay(scRNA.SeuObj_Small) <- "RNA"
#     Idents(seuratObject) <- seuratObject@meta.data[["ReCluster"]]
#     ReCellType.markers <- FindAllMarkers(seuratObject, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
# 
#     write.table(ReCellType.markers, file = paste0( Save.Path,"/ReCelltypeMarker2_AllGene.txt"),
#                 quote = F,sep = "\t",row.names = F)
#     save.image("SeuratObject_CDS_PRJCA001063_MaligAnno.RData")






