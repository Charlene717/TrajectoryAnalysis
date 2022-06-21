## Ref: https://github.com/yanwu2014/swne
## Paper Ref: https://www.sciencedirect.com/science/article/pii/S240547121830440X?via%3Dihub

if(!require(remotes)){ install.packages("remotes")}  # If not already installed; 
  remotes::install_github("linxihui/NNLM")
  remotes::install_github("yanwu2014/swne")
  
  devtools::install_github("aertslab/RcisTarget")
  devtools::install_github("aertslab/AUCell")
  devtools::install_github("aertslab/cisTopic")
  
  
## Gene Expression Quickstart with Seurat
  library(swne)
  ## Load object
  obj <- readRDS("Data/pbmc3k_final.RObj")
  
  ## Extract clusters
  clusters <- obj$seurat_clusters
  
  ## Select genes to embed
  genes.embed <- c("MS4A1", "GNLY", "CD3E", "CD14",
                   "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
  
  ## Run SWNE
  swne.embedding <- RunSWNE(obj, k = 16, genes.embed = genes.embed)
  
  ## Plot SWNE
  PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = clusters,
           do.label = T, label.size = 3.5, pt.size = 1.5, show.legend = F,
           seed = 42)
  
##### PRJCA001063 ##### 
  ## Gene Expression Quickstart with Seurat
  ## Load object
  obj <- seuratObject
  
  ## Extract clusters
  clusters <- seuratObject@meta.data[["Cell_type"]]
  
  ## Select genes to embed
  genes.embed <- c("MS4A1", "GNLY", "CD3E", "CD14",
                   "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
  
  ## Run SWNE
  swne.embedding <- RunSWNE(obj, k = 16, genes.embed = genes.embed)
  
  ## Plot SWNE
  PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = clusters,
           do.label = T, label.size = 3.5, pt.size = 1.5, show.legend = F,
           seed = 42)
    
  