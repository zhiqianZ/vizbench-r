
load_pkgs <- function() {
  library(Seurat)
}

log1pCP10k = function(args){
  message("Running log1pCP10k")
  seurat.obj <- read_seurat(args$simulate.ad)
  NormalizeData(seurat.obj, scale.factor = 10^4)
}

log1pCPM = function(args){
  message("Running log1pCPM")
  seurat.obj <- read_seurat(args$simulate.ad)
  NormalizeData(seurat.obj, scale.factor= 10^6)
}

sctransform = function(args){
  message("Running SCTransform")
  seurat.obj <- read_seurat(args$simulate.ad)
  nhvgs = args$nhvgs
  seurat.obj[["RNA"]] <- split(seurat.obj[["RNA"]], f = seurat.obj$batch)
  seurat.obj = SCTransform(seurat.obj, vst.flavor = "v2", verbose = FALSE, 
              variable.features.n = nhvgs, method = "glmGamPoi",min_cells=0)
  return(seurat.obj)
}












