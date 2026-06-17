
load_pkgs <- function() {
  library(SeuratWrappers)
  library(Seurat)
  library(harmony)
  library(rliger)
  library(Rfast)
  library(dplyr)
  library(anndataR)
  library(rjson)
  library(irlba)
  library(integrateRigor)
}
                
find_hvgs_seuratv5 <- function(seurat.obj, nfeatures = 2000) {
  seurat.obj <- FindVariableFeatures(seurat.obj,
                                     selection.method = "vst", 
                                     nfeatures = nfeatures)
  VariableFeatures(seurat.obj)
}


harmony_integrateRigor = function(args){
  message("Running Harmony with IntegrateRigor")
  nhvgs <- args$nhvgs
  npcs <- args$npcs
  norm_method <- read_normmethod(args$normalize.json)
  obj <- read_seurat(args$normalize.ad)
  obj = BatchStabilityEst(obj, batch="batch", K=5, n.cores=10, subsample=0.03)
  obj = BatchStableGenes(obj, plot = T)
  bsg = obj@misc$batch.stable.genes
  theta = c(2, 4, 8)
  nclust = c(80, 100, 120)
  param = make.parameter.df(theta, nclust)
  obj = obj[bsg, ]
  if(norm_method != "sctransform"){
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)
    hvgs <- find_hvgs_seuratv5(obj, nhvgs)
    VariableFeatures(obj) <- hvgs
    obj <- ScaleData(obj) 
  }else{
    hvgs <- rownames(obj)
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)
    VariableFeatures(obj) <- hvgs
  }
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs)
  obj = IntegrateRigor.ParameterS(obj, method = Harmony, parameter.df = param, force.run = T, ndims.score = npcs, ndims=npcs, subsample=0.03)
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  return(obj)
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
    hvgs <- rownames(so)
    so[["RNA"]] <- split(so[["RNA"]], f = so$batch)
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
    hvgs <- rownames(so)
    so[["RNA"]] <- split(so[["RNA"]], f = so$batch)
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
  #ns <- asNamespace("Seurat")
  #if (bindingIsLocked("RunCCA.default", ns)) unlockBinding("RunCCA.default", ns)
  #assign("RunCCA.default", RunCCA.default2, envir = ns)
  #lockBinding("RunCCA.default", ns)
  #registerS3method("RunCCA", "default", RunCCA.default2, envir = ns)
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
    hvgs <- rownames(so)
    so[["RNA"]] <- split(so[["RNA"]], f = so$batch)
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
    hvgs <- rownames(so)
    so[["RNA"]] <- split(so[["RNA"]], f = so$batch)
    VariableFeatures(so) <- hvgs
    so = IntegrateLayers(
      object = so, method = FastMNNIntegration,
      new.reduction = "integrated", orig.reduction = NULL,
      verbose = TRUE,
      batch = so$batch,
      features = VariableFeatures(so),
      assay.type = "scale.data"
    )
  }
  so <- JoinLayers(so)
  return(so)
}

LIGER = function(args){
  message("Running LIGER")
  nhvgs <- args$nhvgs
  npcs <- args$npcs
  so <- read_seurat(args$simulate.ad)
  
  so <- so %>%
     normalize() %>%
     selectGenes(nGenes = nhvgs) %>%
     scaleNotCenter()
   so <- so %>%
     runINMF(k = npcs) %>%
     quantileNorm()
   so[["integrated"]] = so[["inmfNorm"]]
   so <- JoinLayers(so)
   hvgs = rownames(so[['RNA']]$ligerScaleData)
   so = so[hvgs,]
   return(so)
}


scVI = function(args){
  message("Running scVI")
  nhvgs <- args$nhvgs
  npcs <- args$npcs
  scvi_conda <- args$scvi_conda
  so <- read_seurat(args$simulate.ad)
  
  so[["RNA"]] <- split(so[["RNA"]], f = so$batch)
  hvgs <- find_hvgs_seuratv5(so, nhvgs)
  
  so <- IntegrateLayers(
    so,
    method = scVIIntegration,
    new.reduction = "integrated",
    conda_env = scvi_conda,
    verbose = T,
    features = hvgs,
    ndims = npcs,
    layers = "counts",
    orig.reduction = NULL,
    scale.layer = NULL,
    assay = "RNA"
  )
  so <- JoinLayers(so)
  rownames(so[['integrated']]@cell.embeddings) = colnames(so)
  return(so)
}

