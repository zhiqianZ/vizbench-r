
load_pkgs <- function() {
  library(Seurat)
  library(distances)
  library(dplyr)
}

#use_condaenv("Benchmark")
# MR: comment out python for now; TODO: build/integrate reticulate env
# sc = import("scanpy", convert=F)
# source_python("utils/normalization.py")

# Normalization = function(seurat.obj, NormMethod){
#   seurat.obj = NormMethod(seurat.obj)
#   return(seurat.obj)
# }

# old version
# log1pCP10k = function(seurat.obj){
#   print("Running log1pCP10k")
#   seurat.obj = NormalizeData(seurat.obj,scale.factor = 10^4)
#   return(seurat.obj)
# }

# new version
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
  SCTransform(seurat.obj, vst.flavor = "v2", verbose = FALSE, 
              variable.features.n = nhvgs, method = "glmGamPoi")
}

# log1pPF = function(args){
#   message("Running log1pPF")
#   seurat.obj <- read_seurat(args$simulate.ad)
#   counts.layer = Layers(object = seurat.obj, search = "counts")
#   save <- make.unique(names = gsub(
#     pattern = "counts",
#     replacement = "data",
#     x = counts.layer
#   ))
#   for (i in seq_along(along.with = counts.layer)) {
#     l <- counts.layer[i]
#     LayerData(
#       object = seurat.obj,
#       layer = save[i],
#       features = Features(x = seurat.obj, layer = l),
#       cells = Cells(x = seurat.obj, layer = l)
#     ) <- logPF(t(LayerData(object = seurat.obj, layer = l, fast = NA)))
#   }
#   return(seurat.obj)
# }
# 
# PFlog1pPF = function(args){
#   message("Running PFlog1pPF")
#   seurat.obj <- read_seurat(args$simulate.ad)
#   counts.layer = Layers(object = seurat.obj, search = "counts")
#   save <- make.unique(names = gsub(
#     pattern = "counts",
#     replacement = "data",
#     x = counts.layer
#   ))
#   for (i in seq_along(along.with = counts.layer)) {
#     l <- counts.layer[i]
#     LayerData(
#       object = seurat.obj,
#       layer = save[i],
#       features = Features(x = seurat.obj, layer = l),
#       cells = Cells(x = seurat.obj, layer = l)
#     ) <- PFlogPF(t(LayerData(object = seurat.obj, layer = l, fast = NA)))
#   }
#   return(seurat.obj)
# }

# TODO: rewrite this in native Python
# `log1pCPMedian` = function(seurat.obj){
#   print("Running Scanpy default normalization")
#   counts.layer = Layers(object = seurat.obj, search = "counts")
#   save <- make.unique(names = gsub(
#     pattern = "counts",
#     replacement = "data",
#     x = counts.layer
#   ))
#   for (i in seq_along(along.with = counts.layer)) {
#     l <- counts.layer[i]
#     count = t(LayerData(object = seurat.obj, layer = l, fast = NA))
#     adata = sc$AnnData(count)
#     norm = log1p(py_to_r(sc$pp$normalize_total(adata,inplace=F)$X))
#     LayerData(
#       object = seurat.obj,
#       layer = save[i],
#       features = Features(x = seurat.obj, layer = l),
#       cells = Cells(x = seurat.obj, layer = l)
#     ) <- t(norm)
#   }
#   return(seurat.obj)
# }












