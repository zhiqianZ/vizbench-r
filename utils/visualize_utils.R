load_pkgs <- function() {
  library(Seurat)
}

SeuratUMAP = function(args){
  message("Running SeuratUMAP")
  npcs <- args$npcs
  nthreads <- args$nthreads
  message(args$integrate_raw.ad)
  message(args$integrate_norm.ad)
  if(is.null(args$integrate_norm.ad)){
    fn = args$integrate_raw.ad
  }else{
    fn = args$integrate_norm.ad
  }
  so <- read_seurat(fn)
  so = RunUMAP(so, dims = 1:npcs, 
                       reduction="integrated",
                       n_threads = nthreads)
  Embeddings(so,"umap")
}

BHtSNE = function(args){
  message("Running BH-tSNE")
  message(args$integrate_raw.ad)
  message(args$integrate_norm.ad)
  if(is.null(args$integrate_norm.ad)){
    fn = args$integrate_raw.ad
  }else{
    fn = args$integrate_norm.ad
  }
  so <- read_seurat(fn)
  npcs <- args$npcs
  nthreads <- args$nthreads
  so = RunTSNE(so, dims = 1:npcs, reduction="integrated", 
               tsne.method = "Rtsne",check_duplicates = FALSE,
               num_threads = nthreads)
  Embeddings(so, "tsne")
}

densMAP = function(args){
  message("Running densMAP")
  message(args$integrate_raw.ad)
  message(args$integrate_norm.ad)
  if(is.null(args$integrate_norm.ad)){
    fn = args$integrate_raw.ad
  }else{
    fn = args$integrate_norm.ad
  }
  so <- read_seurat(fn)
  npcs <- args$npcs
  nthreads <- args$nthreads
  so = RunUMAP(so, umap.method = "umap-learn", 
               dims = 1:npcs, reduction="integrated",
               densmap = T, n_jobs = nthreads)
  
  return(Embeddings(so,"umap"))
}

FItSNE =  function(args){
  message("Running FIt-SNE")
  message(args$integrate_raw.ad)
  message(args$integrate_norm.ad)
  if(is.null(args$integrate_norm.ad)){
    fn = args$integrate_raw.ad
  }else{
    fn = args$integrate_norm.ad
  }
  so <- read_seurat(fn)
  npcs <- args$npcs
  nthreads <- args$nthreads
  latent = so@reductions$integrated@cell.embeddings[,1:npcs]
  pcaInit = prcomp(latent, rank=2)$x 
  pcaInit = pcaInit / (sd(pcaInit[,1])* (nrow(pcaInit)-1) / nrow(pcaInit) ) * 0.0001
  vis = fftRtsne(latent, nthreads = nthreads, 
                 perplexity_list = c(30, nrow(latent)/100),  # memory issue, set n/1000 instead of n/100
                 initialization = pcaInit, 
                 learning_rate = nrow(latent)/12)
  rownames(vis) = colnames(so)
  return(vis)
}















