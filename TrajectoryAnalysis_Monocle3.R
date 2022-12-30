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
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "DataSetID", label = F, pt.size = 0.5) 
  
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
             color_cells_by = "Cell_type",  # cluster
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE)
  
  ## Color cells by pseudotime
  # ## Order the cells in pseudotime
  # cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 1]))
  
  ## Order the cells in pseudotime
  cds <- order_cells(cds)
  
  
  plot_cells(cds,
             color_cells_by = "pseudotime",
             group_cells_by = "cluster",
             label_cell_groups = FALSE,
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE,
             label_roots = FALSE,
             trajectory_graph_color = "grey60")
  
  ## Subcluster
  # cds_subset <- choose_cells(cds)  
  cds_sub <- cds[,cds@clusters@listData[["UMAP"]][["clusters"]] %in% c(44)] 
  
  # Subset cells by branch
  cds_sub <- choose_graph_segments(cds)
  cds_sub <- cds[,cds@colData@rownames %in% cds_sub@colData@rownames] 
  
  plot_cells(cds_sub,
             color_cells_by = "cluster",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE)
  
  
  cds_sub <- order_cells(cds_sub)
  plot_cells(cds_sub,
             color_cells_by = "pseudotime",
             group_cells_by = "cluster",
             label_cell_groups = FALSE,
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE,
             label_roots = FALSE,
             trajectory_graph_color = "grey60")
  
#### Insert Monocle3 information to Seurat
  # Error
  # scRNA.SeuObj <- as.Seurat(cds, assay = NULL)
  # FeaturePlot(scRNA.SeuObj, "monocle3_pseudotime")
  
  cds_Meta.df <- data.frame(Cell_ID = names(cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]),
                            Pseudotime = cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]],
                            MonocleCluster =cds@clusters@listData[["UMAP"]][["clusters"]],
                            MonocleClusterP = cds@clusters@listData[["UMAP"]][["partitions"]]) 
  
  cds_SubPseudo.df <- data.frame(Cell_ID = names(cds_sub@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]),
                                 SubPseudotime = cds_sub@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]) 
    
    
  Seurat_Meta.df <- data.frame(Cell_ID = rownames(scRNA.SeuObj@meta.data),
                               scRNA.SeuObj@meta.data %>% as.data.frame())
  
  
  Seurat_Meta.df <- left_join(Seurat_Meta.df, cds_Meta.df, by="Cell_ID")
  Seurat_Meta.df <- left_join(Seurat_Meta.df, cds_SubPseudo.df, by="Cell_ID")
  
  ## Ref: https://github.com/satijalab/seurat/issues/4124
  scRNA.SeuObj@meta.data <- Seurat_Meta.df
  rownames(scRNA.SeuObj@meta.data) <- Seurat_Meta.df$Cell_ID
  
  
  ## UMAP
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "MonocleCluster", label = TRUE, pt.size = 0.5) + NoLegend()
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "MonocleClusterP", label = TRUE, pt.size = 0.5) + NoLegend()
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type", label = TRUE, pt.size = 0.5) + NoLegend()
  
  # scRNA.SeuObj_Ref <- scRNA.SeuObj[,scRNA.SeuObj$DataSetID %in% "PRJCA001063"]
  # memory.limit(700000)
  scRNA_Sub.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj$MonocleClusterP == "1"]
  DimPlot(scRNA_Sub.SeuObj, reduction = "umap",group.by = "Cell_type", label = TRUE, pt.size = 0.5) + NoLegend()
  
  # ## Subset Seurat object based on identity class, also see ?SubsetData
  # ## Ref: https://satijalab.org/seurat/articles/essential_commands.html
  # ## Cholmod error 'out of memory' : Merging Seurat Objects
  # ## Ref: https://stackoverflow.com/questions/66079047/cholmod-error-out-of-memory-merging-seurat-objects
  # 
  # scRNA_Sub.SeuObj <-  subset(x = scRNA.SeuObj, subset = MonocleClusterP == 1)
  # Idents(object = scRNA.SeuObj) <- "Cell_type"
  # subset(x = scRNA.SeuObj, idents = "B_cell")
  
  # ## Error
  # FeaturePlot(scRNA_Sub.SeuObj, reduction = "umap","Pseudotime", label = TRUE, pt.size = 0.5) + NoLegend()+
  #   scale_color_gradient(low = "cyan",high = "red")
  

  
  scRNA_Sub_Bra.SeuObj <- scRNA.SeuObj[,!scRNA.SeuObj$SubPseudotime == Inf]
  DimPlot(scRNA_Sub_Bra.SeuObj, reduction = "umap",group.by = "Cell_type", label = TRUE, pt.size = 0.5) + NoLegend()
  
  scRNA_Sub_Bra.SeuObj <- RunPCA(scRNA_Sub_Bra.SeuObj)
  print(scRNA_Sub_Bra.SeuObj[["pca"]], dims = 1:5, nfeatures = 15)
  
  # scRNA_Sub_Bra.SeuObj[["pca"]]  -> TTT
  
  VizDimLoadings(scRNA_Sub_Bra.SeuObj, dims = 1:2, reduction = "pca")
  DimHeatmap(scRNA_Sub_Bra.SeuObj, dims = 1, cells = 500, balanced = TRUE)
  DimHeatmap(scRNA_Sub_Bra.SeuObj, dims = 1:15, cells = 500, balanced = TRUE)
  AFD_genes <- c("MUC4","SAA4","PPY","MUC21","LYPD2")
  
  
  ## Plot Peudo
  AFD_genes <- c("SAA2", "KRT23", "UBD", "SAA1", "CXCL1", "BIRC3", "CXCL6", "SERPINB2", "LAMC2", "CXCL5")
  AFD_genes <- c("LDHB", "RACK1", "VIM", "MEF2C", "TCF4", "SPARC", "COX7A1", "GMFG", "A2M", "SRGN")

  AFD_genes <- c("SAA2", "KRT23", "UBD", "SOX9", "VNN1", "CX3CL1", "MMP7", "ANXA3", "PDZK1IP1", "TNFRSF12A")
  AFD_genes <- c("MUC4", "MUC21", "LYPD2", "PAEP", "FAM83A", "ADGRF1", "MACROD2", "SAA2", "CXCL5", "KRT81" )
  AFD_genes <- c("CXCL6", "CXCL1", "TFPI2", "CFTR", "SLC4A4", "SERPINA5", "DCDC2", "CCL2", "CXCL8", "SLC34A2","CLU", "SPP1", "VCAM1", "SBSPON", "HSD17B2" )
  
  AFD_genes <- c("CLDN1", "TNFAIP2", "CXCL11", "TACSTD2", "CCL28", "RPL17", "FTL", "BGN", "SPARCL1", "ARHGDIB" ) 
  

  AFD_genes <- c("CLU", "SPP1", "VCAM1", "HSD17B2" )
  AFD_genes <- c("SAA2",  "VNN1","MUC4","ADGRF1", "BGN")
  
  rowData(cds_sub)$gene_short_name <- cds_sub@rowRanges@partitioning@NAMES
  AFD_lineage_cds <- cds_sub[rowData(cds_sub)$gene_short_name %in% AFD_genes,]
  plot_genes_in_pseudotime(AFD_lineage_cds, color_cells_by="Cell_type",cell_size = 2, min_expr=0.1)
  
  
  
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






