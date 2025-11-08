
load_pkgs <- function() {
  library(SeuratWrappers)
  library(Seurat)
  library(harmony)
  #library(rliger)
  library(Rfast)
  library(dplyr)
  library(anndataR)
  library(rjson)
}

find_hvgs_seuratv5 <- function(seurat.obj, nfeatures = 2000) {
  seurat.obj <- FindVariableFeatures(seurat.obj,
                                     selection.method = "vst", 
                                     nfeatures = nfeatures)
  VariableFeatures(seurat.obj)
}

harmony = function(args) {
  message("Running Harmony")
  nhvgs <- args$nhvgs
  npcs <- args$npcs
  norm_method <- read_normmethod(args$normalize.json)
  so <- read_seurat(args$normalize.ad)
  if(norm_method != "sctransform"){
    so[["RNA"]] <- split(so[["RNA"]], f = so$batch)
    hvgs <- find_hvgs_seuratv5(so, nhvgs)
    VariableFeatures(so) <- hvgs
    so <- ScaleData(so) 
  }else{
    hvgs <- read_hvgs(args$sct.hvgs.json)
    VariableFeatures(so) <- hvgs
  }
  so <- RunPCA(so, features = VariableFeatures(so), npcs = npcs)
  so <- RunHarmony(so, "batch", reduction.save = "integrated")
  so <- JoinLayers(so)
  return(so)
}

SeuratRPCA = function(args){
  message("Running Seurat RPCA")
  nhvgs <- args$nhvgs
  npcs <- args$npcs
  norm_method <- read_normmethod(args$normalize.json)
  so <- read_seurat(args$normalize.ad)
  if(norm_method != "sctransform"){
    so[["RNA"]] <- split(so[["RNA"]], f = so$batch)
    hvgs <- find_hvgs_seuratv5(so, nhvgs)
    VariableFeatures(so) <- hvgs
    so <- ScaleData(so)
    so <- RunPCA(so, features = VariableFeatures(so), npcs = npcs)
    so = IntegrateLayers(
       object = so, method = RPCAIntegration,
       orig.reduction = "pca", new.reduction = "integrated",
       verbose = TRUE,
       dims = 1:npcs,
       features = VariableFeatures(so)
     ) 
  }else{
    hvgs <- read_hvgs(args$sct_hvgs.json)
    VariableFeatures(so) <- hvgs
    so <- RunPCA(so, features = VariableFeatures(so), npcs = npcs)
    so = IntegrateLayers(
       object = so, method = RPCAIntegration,
       orig.reduction = "pca", new.reduction = "integrated",
       verbose = TRUE,
       dims = 1:npcs,
       features = VariableFeatures(so)
    )
  }
  so <- JoinLayers(so)
  return(so)
}

SeuratCCA = function(args){
  message("Running Seurat CCA")
  nhvgs <- args$nhvgs
  npcs <- args$npcs
  norm_method <- read_normmethod(args$normalize.json)
  so <- read_seurat(args$normalize.ad)
  if(norm_method != "sctransform"){
    so[["RNA"]] <- split(so[["RNA"]], f = so$batch)
    hvgs <- find_hvgs_seuratv5(so, nhvgs)
    VariableFeatures(so) <- hvgs
    so <- ScaleData(so)
    so <- RunPCA(so, features = VariableFeatures(so), npcs = npcs)
    so = IntegrateLayers(
       object = so, method = CCAIntegration,
       orig.reduction = "pca", new.reduction = "integrated",
       verbose = TRUE,
       dims = 1:npcs,
       features = VariableFeatures(so)
     ) 
  }else{
    hvgs <- read_hvgs(args$sct_hvgs.json)
    VariableFeatures(so) <- hvgs
    so <- RunPCA(so, features = VariableFeatures(so), npcs = npcs)
    so = IntegrateLayers(
       object = so, method = CCAIntegration,
       orig.reduction = "pca", new.reduction = "integrated",
       verbose = TRUE,
       dims = 1:npcs,
       features = VariableFeatures(so)
    )
  }
  so <- JoinLayers(so)
  return(so)
}

fastMNN = function(args) {
  message("Running FastMNN")
  nhvgs <- args$nhvgs
  npcs <- args$npcs
  norm_method <- read_normmethod(args$normalize.json)
  so <- read_seurat(args$normalize.ad)
  if(norm_method != "sctransform"){
    so[["RNA"]] <- split(so[["RNA"]], f = so$batch)
    hvgs <- find_hvgs_seuratv5(so, nhvgs)
    VariableFeatures(so) <- hvgs
    so <- ScaleData(so)
    so <- RunPCA(so, npcs = npcs)
    so = IntegrateLayers(object = so, method = FastMNNIntegration,
                                 new.reduction = 'integrated', verbose = TRUE, 
                                 orig.reduction = NULL,
                                 features = VariableFeatures(so),
                                 assay.type = "logcounts")
  }else{
    hvgs <- read_hvgs(args$sct_hvgs.json)
    VariableFeatures(so) <- hvgs
    so = IntegrateLayers(
      object = so, method = FastMNNIntegration,
      new.reduction = "integrated", orig.reduction = NULL,
      verbose = TRUE,
      batch = so$batch,
      features = VariableFeatures(so),
      assay.type = "scaledata"
    )
  }
  so <- JoinLayers(so)
  return(so)
}

LIGER = function(args){
  message("Running LIGER")
  nhvgs <- args$nhvgs
  npcs <- args$npcs
  norm_method <- read_normmethod(args$normalize.json)
  so <- read_seurat(args$normalize.ad)
  
  so <- so %>%
     normalize() %>%
     selectGenes(nGenes = nhvgs) %>%
     scaleNotCenter()
   so <- so %>%
     runINMF(k = npcs) %>%
     quantileNorm()
   so[["integrated"]] = so[["inmfNorm"]]
   so <- JoinLayers(so)
   return(so)
}


scVI = function(args){
  message("Running scVI")
  message("Running scVI")
  nhvgs <- args$nhvgs
  npcs <- args$npcs
  norm_method <- read_normmethod(args$normalize.json)
  so <- read_seurat(args$normalize.ad)
  
  so[["RNA"]] <- split(so[["RNA"]], f = so$batch)
  hvgs <- find_hvgs_seuratv5(so, nhvgs)
  
  so <- IntegrateLayers(
    so,
    method = scVIIntegration,
    new.reduction = "integrated",
    #conda_env = "Benchmark",
    verbose = T,
    features = hvgs,
    ndims = npcs,
    layers = "counts",
    orig.reduction = NULL,
    scale.layer = NULL,
    assay = "RNA"
  )
  so <- JoinLayers(so)
  return(so)
}

