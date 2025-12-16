
load_pkgs <- function() {
  library(SeuratWrappers)
  library(Seurat)
  library(harmony)
  library(rliger)
  library(Rfast)
  library(dplyr)
  library(anndataR)
  library(rjson)
}

FastMNNIntegration <- function(
    object,
    assay = NULL,
    orig = NULL,
    groups = NULL,
    layers = NULL,
    scale.layer = NULL,
    features = 2000,
    new.reduction = "integrated.mnn",
    reduction.key = "mnn_",
    reconstructed.assay = "mnn.reconstructed",
    verbose = TRUE,
    ...
) {
  print("new fastmnn")
  object <- CreateSeuratObject(object)
  if (is.numeric(x = features)) {
    if (verbose) {
      message(paste("Computing", features, "integration features"))
    }
    features <- SelectIntegrationFeatures5(object = object, features = features)
  }
  if(assay == "SCT"){
    layers <- layers %||% Layers(object, search = "scale.data")  # change to scale.data for sctransform
  }else{
    layers <- layers %||% Layers(object, search = "data")
  }
  if (verbose) {
    message("Converting layers to SingleCellExperiment")
  }
  objects.sce <- lapply(
    X = layers,
    FUN = function(x, f) {
      return(as.SingleCellExperiment(
        x = subset(x = object,
                   features = f,
                   cells = colnames(LayerData(object, layer = x)))
      )
      )
    },
    f = features
  )
  if (verbose) {
    message("Running fastMNN")
  }
  out <- do.call(
    what = batchelor::fastMNN,
    args = c(
      objects.sce,
      list(...)
    )
  )
  colnames(x = SingleCellExperiment::reducedDim(x = out)) <- paste0(reduction.key, 1:ncol(x = SingleCellExperiment::reducedDim(x = out)))
  reduction <- CreateDimReducObject(
    embeddings = SingleCellExperiment::reducedDim(x = out),
    loadings = as.matrix(SingleCellExperiment::rowData(x = out)),
    assay = DefaultAssay(object = object),
    key = reduction.key
  )
  # Add reconstructed matrix (gene x cell)
  reconstructed_assay <- CreateAssayObject(
    data = as.matrix(SummarizedExperiment::assay(x = out)),
  )
  # Add variable features
  VariableFeatures(object = reconstructed_assay) <- features
  #Tool(object = object) <- S4Vectors::metadata(x = out)
  #object <- LogSeuratCommand(object = object)
  output.list <- list(reduction, reconstructed_assay)
  names(output.list) <- c(new.reduction, reconstructed.assay)
  return(output.list)
}


RunCCA.default2 <- function(
  object1,
  object2,
  standardize = TRUE,
  num.cc = 20,
  seed.use = 42,
  verbose = FALSE,
  ...
) {
  print("new")
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  cells1 <- colnames(x = object1)
  cells2 <- colnames(x = object2)
  print("here1")
  message("here1")
  if (standardize) {
    print("here2")
    message("here2")
    object1 <- Standardize(mat = object1, display_progress = FALSE)
    object2 <- Standardize(mat = object2, display_progress = FALSE)
  }
  print("here3")
  message("here3")
  mat3 <- crossprod(x = object1, y = object2)
  print("here4")
  message("here4")
  print(dim(mat3))
  message(dim(mat3))
  cca.svd <- irlba(A = mat3, nv = num.cc)
  print("here5")
  message("here5")
  cca.data <- rbind(cca.svd$u, cca.svd$v)
  colnames(x = cca.data) <- paste0("CC", 1:num.cc)
  rownames(cca.data) <- c(cells1, cells2)
  print("here6")
  message("here6")
  cca.data <- apply(
    X = cca.data,
    MARGIN = 2,
    FUN = function(x) {
      if (sign(x[1]) == -1) {
        x <- x * -1
      }
      return(x)
    }
  )
  print("here7")
  message("here7")
  return(list(ccv = cca.data, d = cca.svd$d))
}

godmode:::assignAnywhere("RunCCA.default", RunCCA.default2)
getFromNamespace("RunCCA.default", "Seurat")
godmode:::assignAnywhere("FastMNNIntegration", FastMNNIntegration)

                
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
  message("Running Seurat CCA")
  getFromNamespace("RunCCA.default", "Seurat")
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

SeuratCCA = function(args){
  message("Running Seurat CCAv4")
  nhvgs <- args$nhvgs
  npcs <- args$npcs
  norm_method <- read_normmethod(args$normalize.json)
  so <- read_seurat(args$normalize.ad)
  if(norm_method != "sctransform"){
    so <- SplitObject(so, split.by = "batch")
    so <- lapply(X = so, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nhvgs)
    })
    anchors <- FindIntegrationAnchors(object.list = so, dims = 1:npcs)
    integrated <- IntegrateData(anchorset = anchors, dims = 1:npcs, new.assay.name = "cca")
    DefaultAssay(integrated) <- "cca"
    integrated <- ScaleData(integrated, verbose = FALSE)
    integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE, reduction.name = "integrated")
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
  return(integrated)
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

