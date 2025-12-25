load_pkgs <- function() {
  library(Seurat)
  library(densvis)
  library(phateR)
}

SeuratUMAP = function(args){
  message("Running SeuratUMAP")
  npcs <- args$npcs
  nthreads <- args$nthreads
  fn = args$integrate.ad
  so <- read_seurat(fn)
  so = RunUMAP(so, dims = 1:npcs, 
                       reduction="integrated",
                       n_threads = nthreads)
  Embeddings(so,"umap")
}

BHtSNE = function(args){
  message("Running BH-tSNE")
  fn = args$integrate.ad
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
  fn = args$integrate.ad
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
  fn = args$integrate.ad
  so <- read_seurat(fn)
  npcs <- args$npcs
  nthreads <- args$nthreads
  latent = Embeddings(so,"integrated")[,1:npcs]
  pcaInit = prcomp(latent, rank=2)$x 
  pcaInit = pcaInit / (sd(pcaInit[,1])* (nrow(pcaInit)-1) / nrow(pcaInit) ) * 0.0001
  vis = fftRtsne(latent, nthreads = nthreads, 
                 perplexity_list = c(30, nrow(latent)/100),  # memory issue, set n/1000 instead of n/100
                 initialization = pcaInit, 
                 learning_rate = nrow(latent)/12)
  rownames(vis) = colnames(so)
  return(vis)
}

denSNE = function(args){
  message("Running denSNE")
  fn = args$integrate.ad
  so <- read_seurat(fn)
  npcs <- args$npcs
  nthreads <- args$nthreads
  
  latent = Embeddings(so,"integrated")[,1:npcs]
  vis = densne(latent, num_threads = nthreads)
  rownames(vis) = colnames(so)
  return(vis)
}

PHATE = function(args){
  message("Running PHATE")
  fn = args$integrate.ad
  so <- read_seurat(fn)
  npcs <- args$npcs
  nthreads <- args$nthreads
  latent = Embeddings(so,"integrated")[,1:npcs]
  vis = phate(latent, n.jobs = nthreads)$embedding
  rownames(vis) = colnames(so)
  return(vis)
}

graphFA = function(args){
  print("Running graphFA")
  ad = import("anndata")
  sc = import("scanpy", convert=F)
  fn = args$integrate.ad
  so <- read_seurat(fn)
  npcs <- args$npcs
  nthreads <- args$nthreads
  
  temp_count = matrix(0, nrow = ncol(so), ncol = nrow(so))
  adata = sc$AnnData(temp_count)
  latent.method.key = "X_integrated"
  adata$obsm[latent.method.key] = Embeddings(so,"integrated")[,1:npcs]
  message("finding neighbors")
  sc$pp$neighbors(adata,use_rep=latent.method.key, n_pcs=as.integer(npcs))
  message("starting layouts")
  sc$tl$draw_graph(adata,layout="fa",n_jobs=as.integer(nthreads))
  message("done")
  vis = as.matrix(adata$obsm$get('X_draw_graph_fa'))
  rownames(vis) = colnames(so)
  return(vis)
}












