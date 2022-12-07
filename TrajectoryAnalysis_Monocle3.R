## Ref: https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/monocle_fixed
## Ref: https://cole-trapnell-lab.github.io/monocle3/docs/differential/

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

  ##### Load Packages #####
  #### Basic and BiocManager installation ####
  source("FUN_Package_InstLoad.R")
  FUN_Basic.set <- c("tidyverse","Seurat","ggplot2","ggpmisc", "stringr","magrittr","dplyr")
  FUN_BiocManager.set <- c("fgsea","AnnotationHub","ensembldb", "SeuratDisk","monocle", "SingleR","scRNAseq","celldex","scran")
  
  FUN_Package_InstLoad(Basic.set = FUN_Basic.set, BiocManager.set = FUN_BiocManager.set)
  rm(FUN_Basic.set, FUN_BiocManager.set)
    
  #### GitHub installation ####
  if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
  library(monocle)
  devtools::install_github("cole-trapnell-lab/garnett")
  library(garnett)
  devtools::install_github('cole-trapnell-lab/monocle3')
  library(monocle3)
  remotes::install_github('satijalab/seurat-wrappers')
  library(SeuratWrappers)
  
  devtools::install_github("thomasp85/patchwork")
  library(patchwork)
   
#### Load data #####
  # load("scRNA.SeuObj_CDS_PRJCA001063.RData")
  load("D:/Dropbox/#_Dataset/Cancer/PDAC/2022-12-05_CTAnno_singleR_RefPRJCA001063_PDAC.RData")

  ## Clean up the object
  rm(list=setdiff(ls(), c("scRNA.SeuObj")))
  
  scRNA.SeuObj$Cell_type <- scRNA.SeuObj$singleR_classic_PredbyscRNA
  scRNA.SeuObj$Cell_type  <- gsub(" ", "_", scRNA.SeuObj$Cell_type)
  
  # ## Sub 
  # scRNA.SeuObj_RMDuc2 <- scRNA.SeuObj[,!scRNA.SeuObj@meta.data[["Cell_type"]] == "Ductal_cell_type_2"]
  
##### Current path and new folder setting* #####
  ProjectName = "TrajAna_PCA"
  Sampletype = "PDAC"

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){ dir.create(Save.Path) }


##### Function setting #####
  ## Call function
  source("FUN_Beautify_ggplot.R") 
    

#### Plot UMAP #####
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type", label = TRUE, pt.size = 0.5) + NoLegend()
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "celltype", label = TRUE, pt.size = 0.5) + NoLegend()
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "ReCluster", label = TRUE, pt.size = 0.5) + NoLegend()
  
  # ## Use Idents to plot 
  #   Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["Cell_type"]]
  #   DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  # 
  #   Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["DataSetID"]]
  #   DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  #   FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))
  

## Ref: https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/monocle_fixed
##  Setting up monocle3 cell_data_set object using the SueratWrappers
  library(SeuratWrappers)
  cds <- as.cell_data_set(scRNA.SeuObj)
  cds <- cluster_cells(cds, resolution=1e-3)
  
  p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
  p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
  # p1+p2
  
  library(patchwork)
  wrap_plots(p1, p2)
  
  
  ## Trajectory analysis
  cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
  plot_cells(cds,
             color_cells_by = "cluster",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE)
  
  
#### Save RData #####    
  save.image(paste0(Save.Path,"/scRNA.SeuObj_CDS_PRJCA001063_TrajAna_PCA.RData"))
    
# #### Cell type markers #####
#     ## Create cell type markers dataframe
#     # DefaultAssay(scRNA.SeuObj_Small) <- "RNA"
#     Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["ReCluster"]]
#     ReCellType.markers <- FindAllMarkers(scRNA.SeuObj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
# 
#     write.table(ReCellType.markers, file = paste0( Save.Path,"/ReCelltypeMarker2_AllGene.txt"),
#                 quote = F,sep = "\t",row.names = F)
#     save.image("scRNA.SeuObj_CDS_PRJCA001063_MaligAnno.RData")






