load_pkgs <- function() {
  library(Seurat)
  library(densvis)
  library(phateR)
  library(scDEED)
  library(pracma)
}

## ---- Global ----------------------------------------------------------
SUBSAMPLE_N <- 30000   # cells used for hyperparameter selection (all methods)
SEED        <- 42

## ---- Hyperparameter grids --------------------------------------------------
## min.dist: 3 values, as requested.
GRID_UMAP <- expand.grid(
  n.neighbors = c(5, 20, 30, 40, 50),
  min.dist    = c(0.1, 0.3, 0.5),
  stringsAsFactors = FALSE
)
GRID_SCANPY_UMAP <- expand.grid(
  n_neighbors = c(5, 20, 30, 40, 50),
  min_dist    = c(0.1, 0.3, 0.5),
  stringsAsFactors = FALSE
)

GRID_TSNE <- data.frame(
  perplexity = seq(from = 20, to = 410, by = 30)
)

## densMAP: UMAP grid crossed with a few dens_lambda values.
## dens_lambda controls how strongly local density is preserved 
GRID_DENSMAP <- expand.grid(
  n.neighbors = c(5, 20, 30, 40, 50),
  min.dist    = c(0.1, 0.3, 0.5),
  dens_lambda = c(0.5, 2.0),
  stringsAsFactors = FALSE
)

GRID_GRAPHFA <- data.frame(
  n_neighbors = c(5, 15, 30, 50)
)

SeuratUMAP <- function(args) {
  message("Running SeuratUMAP")
  so <- read_seurat(args$integrate.ad)
  so <- RunUMAP(so, dims = 1:args$npcs, reduction = "integrated",
                umap.method = "umap-learn", n_jobs = args$nthreads)
  Embeddings(so, "umap")
}

BHtSNE <- function(args) {
  message("Running BH-tSNE")
  so <- read_seurat(args$integrate.ad)
  so <- RunTSNE(so, dims = 1:args$npcs, reduction = "integrated",
                tsne.method = "Rtsne", check_duplicates = FALSE,
                num_threads = args$nthreads)
  Embeddings(so, "tsne")
}

densMAP <- function(args) {
  message("Running densMAP")
  so <- read_seurat(args$integrate.ad)
  so <- RunUMAP(so, dims = 1:args$npcs, reduction = "integrated",
                umap.method = "umap-learn", densmap = TRUE,
                n_jobs = args$nthreads)
  Embeddings(so, "umap")
}

FItSNE <- function(args) {
  message("Running FIt-SNE")
  so     <- read_seurat(args$integrate.ad)
  latent <- Embeddings(so, "integrated")[, 1:args$npcs]

  pcaInit <- prcomp(latent, rank = 2)$x
  pcaInit <- pcaInit / (sd(pcaInit[, 1]) * (nrow(pcaInit) - 1) / nrow(pcaInit)) * 1e-4

  vis <- fftRtsne(
    latent,
    nthreads         = args$nthreads,
    perplexity_list  = c(30, nrow(latent) / 100),  # multi-scale; NOT optimized
    initialization   = pcaInit,
    learning_rate    = nrow(latent) / 12
  )
  rownames(vis) <- colnames(so)
  vis
}

denSNE <- function(args) {
  message("Running denSNE")
  so     <- read_seurat(args$integrate.ad)
  latent <- Embeddings(so, "integrated")[, 1:args$npcs]
  vis <- densne(latent, num_threads = args$nthreads)
  rownames(vis) <- colnames(so)
  vis
}

PHATE <- function(args) {
  message("Running PHATE")
  so     <- read_seurat(args$integrate.ad)
  latent <- Embeddings(so, "integrated")[, 1:args$npcs]
  vis <- phate(latent, n.jobs = args$nthreads)$embedding
  rownames(vis) <- colnames(so)
  vis
}

.graphfa_embed <- function(so, npcs, nthreads, n_neighbors) {
  sc <- reticulate::import("scanpy", convert = FALSE)

  temp_count <- matrix(0, nrow = ncol(so), ncol = nrow(so))
  adata      <- sc$AnnData(temp_count)
  adata$obsm["X_integrated"] <- Embeddings(so, "integrated")[, 1:npcs]

  sc$pp$neighbors(adata, use_rep = "X_integrated",
                  n_pcs       = as.integer(npcs),
                  n_neighbors = as.integer(n_neighbors))
  sc$tl$draw_graph(adata, layout = "fa", n_jobs = as.integer(nthreads))

  vis <- as.matrix(adata$obsm$get("X_draw_graph_fa"))
  rownames(vis) <- colnames(so)
  vis
}

graphFA <- function(args) {
  message("Running graphFA")
  so <- read_seurat(args$integrate.ad)
  vis <- .graphfa_embed(so, npcs = args$npcs, nthreads = args$nthreads,
                        n_neighbors = 15)   # scanpy default
  vis
}

permute_reduction <- function(so, reduction = "integrated", key = "integrated_",
                              seed = SEED) {
  emb  <- Embeddings(so, reduction)
  perm <- emb
  set.seed(seed)                          
  for (i in seq_len(ncol(emb))) {
    perm[, i] <- emb[pracma::randperm(nrow(emb)), i]
  }
  so.permuted <- so
  so.permuted[[reduction]] <- CreateDimReducObject(
    embeddings = perm,
    key        = key,
    assay      = DefaultAssay(so)
  )
  so.permuted
}

scdeed_optimize <- function(so, so.permuted, grid, embed_fn,
                            npcs, pre_embedding = "integrated",
                            similarity_percent = 0.5,
                            dubious_cutoff = 0.05, trustworthy_cutoff = 0.95) {

  ## Pre-embedding (high-dimensional) distances -- computed ONCE, reused for
  ## every hyperparameter setting. This is the expensive part of scDEED and
  ## does not depend on the embedder.
  pre <- scDEED::Distances.pre_embedding(
    so, so.permuted, K = npcs, pre_embedding = pre_embedding
  )

  n_dubious <- integer(nrow(grid))

  for (i in seq_len(nrow(grid))) {
    params <- as.list(grid[i, , drop = FALSE])
    message(sprintf("  [%d/%d] %s", i, nrow(grid),
                    paste(names(params), unlist(params), sep = "=", collapse = ", ")))

    d_orig <- distances::distances(embed_fn(so,          params))
    d_perm <- distances::distances(embed_fn(so.permuted, params))

    sim <- scDEED::Cell.Similarity(
      pre$pre_embedding_distances,
      pre$pre_embedding_distances_permuted,
      d_orig, d_perm,
      similarity_percent = similarity_percent
    )
    cls <- scDEED::Cell.Classify(
      sim$rho_original, sim$rho_permuted,
      dubious_cutoff     = dubious_cutoff,
      trustworthy_cutoff = trustworthy_cutoff
    )
    n_dubious[i] <- length(cls$dubious_cells)
  }

  out  <- cbind(grid, number_dubious_cells = n_dubious)
  best <- out[which.min(out$number_dubious_cells), , drop = FALSE]

  message("Best setting:")
  print(best)
  list(num_dubious = out, best = best)
}

prep_for_scdeed <- function(so, n = SUBSAMPLE_N, reduction = "integrated",
                            seed = SEED) {
  n <- min(n, ncol(so))                   # don't oversample small objects
  set.seed(seed)
  sub <- so[, sample(seq_len(ncol(so)), n)]
  list(sub = sub, permuted = permute_reduction(sub, reduction, seed = seed))
}

.embed_umap <- function(so, params, npcs, nthreads) {
  so <- RunUMAP(so, dims = 1:npcs, reduction = "integrated",
                umap.method  = "umap-learn",
                n.neighbors  = params$n.neighbors,
                min.dist     = params$min.dist,
                n_jobs       = nthreads,
                verbose = FALSE)
  Embeddings(so, "umap")
}
.embed_graphfa <- function(so, params, npcs, nthreads) {
  .graphfa_embed(so, npcs = npcs, nthreads = nthreads,
                 n_neighbors = params$n_neighbors)
}
.embed_tsne <- function(so, params, npcs, nthreads) {
  so <- RunTSNE(so, dims = 1:npcs, reduction = "integrated",
                tsne.method      = "Rtsne",
                check_duplicates = FALSE,      # matches the BHtSNE baseline
                num_threads      = nthreads,
                perplexity       = params$perplexity)
  Embeddings(so, "tsne")
}
.embed_densmap <- function(so, params, npcs, nthreads) {
  so <- RunUMAP(so, dims = 1:npcs, reduction = "integrated",
                umap.method  = "umap-learn",
                densmap      = TRUE,
                dens_lambda  = params$dens_lambda,
                n.neighbors  = params$n.neighbors,
                min.dist     = params$min.dist,
                n_jobs       = nthreads,
                verbose = FALSE)
  Embeddings(so, "umap")
}
.embed_scanpy_umap <- function(so, params, npcs, nthreads) {
  latent <- Embeddings(so, "integrated")[, 1:npcs, drop = FALSE]

  vis <- scanpy_umap_from_matrix(
    latent      = latent,
    n_neighbors = as.integer(params$n_neighbors),
    min_dist    = as.numeric(params$min_dist),
    n_jobs      = as.integer(nthreads),
    seed        = 100L
  )
  vis <- as.matrix(vis)
  rownames(vis) <- colnames(so)
  vis
}


SeuratUMAP_scDEED <- function(args) {
  message("Running SeuratUMAP + scDEED")
  npcs     <- args$npcs
  nthreads <- args$nthreads
  so_full  <- read_seurat(args$integrate.ad)

  p   <- prep_for_scdeed(so_full)
  res <- scdeed_optimize(
    p$sub, p$permuted, GRID_UMAP,
    embed_fn = function(so, params) .embed_umap(so, params, npcs, nthreads),
    npcs = npcs
  )
  rm(p); gc()

  best <- res$best
  so_full <- RunUMAP(so_full, dims = 1:npcs, reduction = "integrated",
                     umap.method = "umap-learn",
                     n.neighbors = best$n.neighbors,
                     min.dist    = best$min.dist,
                     n_jobs      = nthreads)

  structure(Embeddings(so_full, "umap"), scdeed = res$num_dubious)
}


BHtSNE_scDEED <- function(args) {
  message("Running BH-tSNE + scDEED")
  npcs     <- args$npcs
  nthreads <- args$nthreads
  so_full  <- read_seurat(args$integrate.ad)

  p   <- prep_for_scdeed(so_full)
  res <- scdeed_optimize(
    p$sub, p$permuted, GRID_TSNE,
    embed_fn = function(so, params) .embed_tsne(so, params, npcs, nthreads),
    npcs = npcs
  )
  rm(p); gc()

  best <- res$best
  so_full <- RunTSNE(so_full, dims = 1:npcs, reduction = "integrated",
                     tsne.method      = "Rtsne",
                     check_duplicates = FALSE,
                     num_threads      = nthreads,
                     perplexity       = best$perplexity)

  structure(Embeddings(so_full, "tsne"), scdeed = res$num_dubious)
}


densMAP_scDEED <- function(args) {
  message("Running densMAP + scDEED")
  npcs     <- args$npcs
  nthreads <- args$nthreads
  so_full  <- read_seurat(args$integrate.ad)

  p   <- prep_for_scdeed(so_full)
  res <- scdeed_optimize(
    p$sub, p$permuted, GRID_DENSMAP,
    embed_fn = function(so, params) .embed_densmap(so, params, npcs, nthreads),
    npcs = npcs
  )
  rm(p); gc()

  best <- res$best
  so_full <- RunUMAP(so_full, dims = 1:npcs, reduction = "integrated",
                     umap.method = "umap-learn",
                     densmap     = TRUE,
                     dens_lambda = best$dens_lambda,
                     n.neighbors = best$n.neighbors,
                     min.dist    = best$min.dist,
                     n_jobs      = nthreads)

  structure(Embeddings(so_full, "umap"), scdeed = res$num_dubious)
}


graphFA_scDEED <- function(args) {
  message("Running graphFA + scDEED")
  npcs     <- args$npcs
  nthreads <- args$nthreads
  so_full  <- read_seurat(args$integrate.ad)

  p   <- prep_for_scdeed(so_full)
  res <- scdeed_optimize(
    p$sub, p$permuted, GRID_GRAPHFA,
    embed_fn = function(so, params) .embed_graphfa(so, params, npcs, nthreads),
    npcs = npcs
  )
  rm(p); gc()

  best <- res$best
  vis  <- .graphfa_embed(so_full, npcs = npcs, nthreads = nthreads,
                         n_neighbors = best$n_neighbors)

  structure(vis, scdeed = res$num_dubious)
}

scanpyUMAP_scDEED <- function(args) {
  message("Running scanpyUMAP + scDEED")
  npcs     <- args$npcs
  nthreads <- args$nthreads
  so_full  <- read_seurat(args$integrate.ad)

  p   <- prep_for_scdeed(so_full)
  res <- scdeed_optimize(
    p$sub, p$permuted, GRID_SCANPY_UMAP,
    embed_fn = function(so, params) .embed_scanpy_umap(so, params, npcs, nthreads),
    npcs = npcs
  )
  rm(p); gc()

  best   <- res$best
  latent <- Embeddings(so_full, "integrated")[, 1:npcs, drop = FALSE]

  vis <- scanpy_umap_from_matrix(
    latent      = latent,
    n_neighbors = as.integer(best$n_neighbors),
    min_dist    = as.numeric(best$min_dist),
    n_jobs      = as.integer(nthreads),
    seed        = 42L
  )
  vis <- as.matrix(vis)
  rownames(vis) <- colnames(so_full)

  structure(vis, scdeed = res$num_dubious)
}


#Distances.UMAP_mod = function(pbmc,pbmc.permuted, K, pre_embedding = 'pca', n = 30, m = 0.3, rerun = T) {
#  distances <- distances::distances
#  if(rerun){
#    pbmc <- Seurat::RunUMAP(pbmc, dims = 1:K, seed.use = 100, reduction = pre_embedding, n.neighbors = n, min.dist = m,
#                       umap.method = "umap-learn")
#}
  
#  UMAP_distances = distances(pbmc@reductions$umap@cell.embeddings)
#  pbmc.permuted <- Seurat::RunUMAP(pbmc.permuted, dims = 1:K, seed.use = 100, 
#                                   reduction = pre_embedding, n.neighbors = n, min.dist = m, 
#                                   umap.method = "umap-learn")
#  UMAP_distances_permuted = distances(pbmc.permuted@reductions$umap@cell.embeddings)
#  results.PCA   <- list("reduced_dim_distances" = UMAP_distances, "reduced_dim_distances_permuted" = UMAP_distances_permuted)
  
#  return(results.PCA)
#}

#assignInNamespace(
#  x = "Distances.UMAP",
#  value = Distances.UMAP_mod,
#  ns = "scDEED"
#)













