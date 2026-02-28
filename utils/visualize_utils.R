load_pkgs <- function() {
  library(Seurat)
  library(densvis)
  library(phateR)
  library(scDEED)
  library(pracma)
}

SeuratUMAP = function(args){
  message("Running SeuratUMAP")
  npcs <- args$npcs
  nthreads <- args$nthreads
  fn = args$integrate.ad
  so <- read_seurat(fn)
  so = RunUMAP(so, dims = 1:npcs, 
                   reduction="integrated",
                   umap.method = "umap-learn")
  Embeddings(so, "umap")
}

BHtSNE = function(args){
  message("Running BH-tSNE")
  fn = args$integrate.ad
  so <- read_seurat(fn)
  npcs <- args$npcs
  nthreads <- args$nthreads
  so = RunTSNE(so, dims = 1:npcs, reduction="integrated", 
               tsne.method = "Rtsne", check_duplicates = FALSE,
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
  
  return(Embeddings(so, "umap"))
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
  sc$pp$neighbors(adata,use_rep=latent.method.key, n_pcs=as.integer(npcs))
  sc$tl$draw_graph(adata,layout="fa",n_jobs=as.integer(nthreads))
  vis = as.matrix(adata$obsm$get('X_draw_graph_fa'))
  rownames(vis) = colnames(so)
  return(vis)
}

SeuratUMAP_scDEED = function(args){
  message("Running SeuratUMAP")
  npcs <- args$npcs
  nthreads <- args$nthreads
  fn = args$integrate.ad
  so_origin <- read_seurat(fn)
  set.seed(1000)
  so = so_origin[, sample(1:ncol(so_origin), 30000)]

  pre_embedding = 'integrated'
  integrated_space = Embeddings(so, pre_embedding)
  permuted_space = integrated_space


  for (i in 1:dim(integrated_space)[2]){
    row = randperm(dim(permuted_space)[1])
    permuted_space[,i]=integrated_space[row,i]
  }
  so.permuted = so
  so.permuted[[pre_embedding]] <- CreateDimReducObject(embeddings = permuted_space , key = "integrated_", assay = DefaultAssay(so.permuted))
  start = Sys.time()
  result = scDEED(so, K = npcs, pre_embedding = pre_embedding, permuted = so.permuted, reduction.method = 'umap')
  for (i in 1:dim(integrated_space)[2]){
    row = randperm(dim(permuted_space)[1])
    permuted_space[,i]=integrated_space[row,i]
  }
  so.permuted = so
  so.permuted[[pre_embedding]] <- CreateDimReducObject(embeddings = permuted_space , key = "integrated_", assay = DefaultAssay(so.permuted))

  result = scDEED(so, K = npcs, pre_embedding = pre_embedding, permuted = so.permuted, reduction.method = 'umap')

  n.neighbors = result$num_dubious[which.min(result$num_dubious$number_dubious_cells), "n_neighbors"]
  min.dist = result$num_dubious[which.min(result$num_dubious$number_dubious_cells), "min.dist"]
  rm(so, so.permuted)
  so_origin = RunUMAP(so_origin, dims = 1:50, reduction="integrated", umap.method = "umap-learn",
               min.dist = min.dist, n.neighbors = n.neighbors)
  
  Embeddings(so_origin, "umap")
}

BHtSNE_scDEED = function(args){
  message("Running SeuratUMAP")
  npcs <- args$npcs
  nthreads <- args$nthreads
  fn = args$integrate.ad
  so_origin <- read_seurat(fn)
  set.seed(1000)
  so = so_origin[, sample(1:ncol(so_origin), 10000)]

  pre_embedding = 'integrated'
  integrated_space = Embeddings(so, pre_embedding)
  permuted_space = integrated_space
  set.seed(1000)
  for (i in 1:dim(integrated_space)[2]){
    row = randperm(dim(permuted_space)[1])
    permuted_space[,i]=integrated_space[row,i]
  }
  so.permuted = so
  so.permuted[[pre_embedding]] <- CreateDimReducObject(embeddings = permuted_space , key = "integrated_", assay = DefaultAssay(so.permuted))

  result = scDEED(so, K = npcs, pre_embedding = pre_embedding, permuted = so.permuted, reduction.method = 'tsne')

  perlexity = result$num_dubious[which.min(result$num_dubious$number_dubious_cells), "perplexity"]
  rm(so, so.permuted)
  so_origin = RunTSNE(so_origin, dims = 1:npcs, reduction="integrated", 
               perplexity = perplexity,
               tsne.method = "Rtsne", check_duplicates = FALSE,
               num_threads = nthreads)
  Embeddings(so_origin, "tsne")
}

Distances.UMAP_mod = function(pbmc,pbmc.permuted, K, pre_embedding = 'pca', n = 30, m = 0.3, rerun = T) {
  distances <- distances::distances
  if(rerun){
    pbmc <- Seurat::RunUMAP(pbmc, dims = 1:K, seed.use = 100, reduction = pre_embedding, n.neighbors = n, min.dist = m,
                       umap.method = "umap-learn")
}
  
  UMAP_distances = distances(pbmc@reductions$umap@cell.embeddings)
  pbmc.permuted <- Seurat::RunUMAP(pbmc.permuted, dims = 1:K, seed.use = 100, 
                                   reduction = pre_embedding, n.neighbors = n, min.dist = m, 
                                   umap.method = "umap-learn")
  UMAP_distances_permuted = distances(pbmc.permuted@reductions$umap@cell.embeddings)
  results.PCA   <- list("reduced_dim_distances" = UMAP_distances, "reduced_dim_distances_permuted" = UMAP_distances_permuted)
  
  return(results.PCA)
}

assignInNamespace(
  x = "Distances.UMAP",
  value = Distances.UMAP_mod,
  ns = "scDEED"
)













