
load_pkgs <- function() {
  library(Seurat)
  library(SanityR)
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
              variable.features.n = nhvgs, method = "glmGamPoi")
  
  seurat.obj = seurat.obj[rownames(seurat.obj[['SCT']]$scale.data),]
  return(seurat.obj)
}


Sanity = function(args){
  seurat.obj <- read_seurat(args$simulate.ad)
  sce <- read_sce(args$simulate.ad)
  # run sanity by batch 
  sce_list <- split(seq_len(ncol(sce)), sce$batch)
  sce_list <- lapply(sce_list, function(idx) {
  sub <- sce[, idx]
  sf <- scater::librarySizeFactors(sub)
  sizeFactors(sub) <- sf / mean(sf)
  sub <- Sanity(sub)
  sub
  })
  norm_mat <- logcounts(sce_sanity)
  # Ensure cell order matches Seurat object
  norm_mat <- norm_mat[, colnames(seurat.obj)]
  
  seurat.obj <- SetAssayData(
  object = seurat.obj,
  assay = "RNA",
  slot = "data",
  new.data = norm_mat
  )
  seurat.obj
}









