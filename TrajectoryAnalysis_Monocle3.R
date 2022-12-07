## Ref: https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/monocle_fixed
## Ref: https://cole-trapnell-lab.github.io/monocle3/docs/differential/

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
  
  
  

    
#### Save RData #####    
  save.image(paste0(Save.Path,"/scRNA.SeuObj_CDS_PRJCA001063_TrajAna.RData"))
    
# #### Cell type markers #####
#     ## Create cell type markers dataframe
#     # DefaultAssay(scRNA.SeuObj_Small) <- "RNA"
#     Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["ReCluster"]]
#     ReCellType.markers <- FindAllMarkers(scRNA.SeuObj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
# 
#     write.table(ReCellType.markers, file = paste0( Save.Path,"/ReCelltypeMarker2_AllGene.txt"),
#                 quote = F,sep = "\t",row.names = F)
#     save.image("scRNA.SeuObj_CDS_PRJCA001063_MaligAnno.RData")






