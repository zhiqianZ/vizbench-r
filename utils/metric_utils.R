
load_pkgs <- function() {
  library(SeuratWrappers)
  library(Seurat)
  library(irlba)
  library(readr)
  library(parallel)
  library(lisi)
# library(distances)
# library(cluster)
# library(Rfast)
# library(dplyr)
}

  
celltype_shape = function(args) {
  
  # read embeddings
  data <- read_csv(args$visualize.csv.gz)
  # read SCE to get batch/celltype 
  sce <- read_sce(args$integrate.ad)
  batch <- sce$batch
  celltype <- sce$celltype
  rm(sce)
  
  val = sapply(unique(batch), function(b){
    tapply(which(batch==b), celltype[batch==b], function(id){
      if(ncol(data)>2){
        svd_res = irlba(apply(data[id,], 2, function(d){d-mean(d)}), nv=2)
      }else{
        svd_res = svd(apply(data[id,], 2, function(d){d-mean(d)}))
      }
      as.numeric((svd_res$d[2] / svd_res$d[1]) >= 0.25)
    })
  })
  return(mean(val,na.rm=TRUE))
}


batch_mixture <- function(args, seed=42, B = 100, n = 10000){
  nthreads = args$nthreads
  set.seed(seed)
  
  # read embeddings
  data <- read_csv(args$visualize.csv.gz)
  # read SCE to get batch/celltype 
  sce <- read_sce(args$integrate.ad)
  batch <- sce$batch
  celltype <- sce$celltype
  rm(sce)
  
  if(nrow(data)<100000){
    B = 1
    n = nrow(data)
    nthreads = 1
  }
  
  val = mclapply(1:B, FUN=function(i){
    id = sample(1:nrow(data), n, replace = FALSE)
    res = mean(sapply(unique(celltype[id]), function(c){
      ids = id[celltype[id]==c]
      if(length(unique(batch[ids]))>1){
        per = min(round(table(batch[ids])/2),30)
        lisi_res  = compute_lisi(data[ids,], data.frame(batch=batch[ids]), 
                                 "batch", perplexity = per)
        return(mean(lisi_res$batch))
      }else{
        return(NA)
      }
    }
    ),na.rm=TRUE)
    return(res)
  }, mc.cores = nthreads)
  val = unlist(val)

  return(mean(val,na.rm=T))
}


celltype_separation = function(args, seed=42 B = 100, n = 10000){
  nthreads = args$nthreads
  set.seed(seed)
  
  # read embeddings
  data <- read_csv(args$visualize.csv.gz)
  # read SCE to get batch/celltype 
  sce <- read_sce(args$integrate.ad)
  batch <- sce$batch
  celltype <- sce$celltype
  rm(sce)
 
  if(nrow(data)<100000){
    nthreads = 1
    B = 1
    n = nrow(data)
  }
  
  val = mclapply(1:B, FUN=function(i){
    id = sample(1:nrow(data), n, replace = F)
    res = mean(sapply(unique(batch[id]), function(b){
      ids = id[batch[id]==b]
      if(length(unique(celltype[ids]))>1){
        dist = distances(data[ids,])
        sil_res = silhouette(as.numeric(celltype[ids]), dist)
        return(mean(tapply(sil_res[,"sil_width"], celltype[ids], mean)))
      }else{
        return(NA)
      }
    }
    ), na.rm=T)
    return(res)
  },mc.cores = n.cores)
  val = unlist(val)
  return(mean(val,na.rm=T))
}



DistancePreservation = function(data, mean_par, celltype, batch, seed=42, dist = "pca"){
  val = sapply(unique(batch), function(b){
    h = mean_par[[b]]
    l = apply(data[batch==b,],2, function(d) tapply(d, celltype[batch==b], mean))
    h = h[complete.cases(h),]
    l = l[complete.cases(l),]
    if(dist == "pca"){
      v = apply(h, 2, sd)
      h = h[,v!=0]
      pca_res <- prcomp(h, center = TRUE, scale. = TRUE)
      pca_scores <- pca_res$x
      dh = as.vector(distances::distance_matrix(distances::distances(pca_scores)))
    }else{
      dh = as.vector(distances::distance_matrix(distances::distances(h)))
    }
    if(dist == "both-pca"){
      v = apply(l, 2, sd)
      l = l[,v!=0]
      pca_res <- prcomp(l, center = TRUE, scale. = TRUE)
      pca_scores <- pca_res$x
      dl = as.vector(distances::distance_matrix(distances::distances(pca_scores)))
    }else{
      dl = as.vector(distances::distance_matrix(distances::distances(l))) 
    }
    cor(dh, dl, method = "spearman")
  })
  return(mean(val))
}

VariancePreservation = function(data, variance_par, celltype, batch, seed=42, dist = "raw"){
  val = sapply(unique(batch), function(b){
    
    h = mean_par[[b]]
    l = apply(data[batch==b,],2, function(d) tapply(d, celltype[batch==b], var))
    h = h[complete.cases(h),]
    l = l[complete.cases(l),]
    if(dist=="pca"){
      v = apply(h, 2, sd)
      h = h[,v!=0]
      pca_res <- prcomp(h, center = TRUE, scale. = TRUE)
      pca_scores <- pca_res$x
      hv = apply(pca_scores, 1, var)
    }else{
      hv = rowSums(h)
    }
    if(dist=="both-pca"){
      v = apply(l, 2, sd)
      l = l[,v!=0]
      pca_res <- prcomp(l, center = TRUE, scale. = TRUE)
      pca_scores <- pca_res$x
      lv = apply(pca_scores, 1, var)
    }else{
      lv = rowSums(l)
    }
    cor(hv, lv, method = "spearman")
  })
  return(mean(val))
}

VarianceSampleSize = function(data, variance_par, celltype, batch, seed=42){
  # h = tapply(batch,celltype,length)
  # l = rowSums(apply(data,2, function(d) tapply(d, celltype, var)))
  # h = h[complete.cases(h)]
  # l = l[complete.cases(l)]
  # val = 1 - abs(cor(h, l, method = "spearman"))
  val = sapply(unique(batch), function(b){
    h = tapply(which(batch==b), celltype[batch==b], length)
    l = rowSums(apply(data[batch==b,],2, function(d) tapply(d, celltype[batch==b], var)))
    h = h[complete.cases(h)]
    l = l[complete.cases(l)]
    1 - abs(cor(h, l, method = "spearman"))
  })
  return(mean(val))
}

LibrarySize = function(data, library_size, celltype, batch, seed=42, B, n, n.cores){
  set.seed(seed)
  if(B == 1){
    n.cores=1
  }
  val = mclapply(1:B, FUN=function(i){
    id = sample(1:nrow(data), n, replace = F)
    res = mean(sapply(unique(batch[id]), function(b){
      ids = id[batch[id]==b]
      if(length(unique(celltype[ids]))>1){
        return(1 - tapply(ids, celltype[ids], function(idd) dcor(library_size[idd], data[idd,])$dcor))
      }else{
        return(NA)
      }
    }),na.rm=T)
    return(res)
  },mc.cores = n.cores)
  val = unlist(val)
  
  return(mean(val))
}

ZeroProportion = function(data, zero_proportion, celltype, batch, seed=42, B, n, n.cores){
  set.seed(seed)
  if(B == 1){
    n.cores=1
  }
  val = mclapply(1:B, FUN=function(i){
    id = sample(1:nrow(data), n, replace = F)
    res = mean(sapply(unique(batch[id]), function(b){
      ids = id[batch[id]==b]
      if(length(unique(celltype[ids]))>1){
        return(1 - tapply(ids, celltype[ids], function(idd) Rfast::dcor(zero_proportion[idd], data[idd,])$dcor))
      }else{
        return(NA)
      }
    }),na.rm=T)
    return(res)
  },mc.cores = n.cores)
  val = unlist(val)

  return(mean(val))
}
