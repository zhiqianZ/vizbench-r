
load_pkgs <- function() {
  library(SeuratWrappers)
  library(Seurat)
  library(SingleCellExperiment)
  library(irlba)
  library(readr)
  library(parallel)
  library(lisi)
  library(distances)
  library(cluster)
  library(Rfast)
  library(dplyr)
}

  
celltype_shape = function(args) {
  
  # read embeddings
  data <- as.data.frame(read_csv(args$visualize.csv.gz))
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


batch_mixture <- function(args, seed=42){
  nthreads = args$nthreads
  B = args$B
  n = args$N
  set.seed(seed)
  
  # read embeddings
  data <- as.data.frame(read_csv(args$visualize.csv.gz))
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


celltype_separation = function(args, seed=42){
  nthreads = args$nthreads
  B = args$B
  n = args$N
  set.seed(seed)
  
  # read embeddings
  data <- as.data.frame(read_csv(args$visualize.csv.gz))
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
  }, mc.cores = nthreads)
  val = unlist(val)
  return(mean(val,na.rm=T))
}



distance_preservation = function(args, seed=42){
  nthreads = args$nthreads
  set.seed(seed)
  
  # read embeddings
  data <- as.data.frame(read_csv(args$visualize.csv.gz))
  mean_par <- as.data.frame(read_csv(args$simulate_mean.csv.gz))
  # read SCE to get batch/celltype 
  sce <- read_sce(args$integrate.ad)
  batch <- sce$batch
  celltype <- sce$celltype
  rm(sce)
  
  val = sapply(unique(batch), function(b){
    h = mean_par[mean_par$batch==b,]
    l = apply(data[batch==b,],2, function(d) tapply(d, celltype[batch==b], mean))
    l = l[complete.cases(l),]
    l = l[h$celltype,]
    h = h[complete.cases(h), 3:ncol(h)]
    v = apply(h, 2, sd)
    h = h[,v!=0]
    pca_res <- prcomp(h, center = TRUE, scale. = TRUE)
    pca_scores <- pca_res$x
    dh = as.vector(distances::distance_matrix(distances::distances(pca_scores)))
    
    dl = as.vector(distances::distance_matrix(distances::distances(l))) 
    
    cor(dh, dl, method = "spearman")
  })
  return(mean(val))
}

variance_preservation = function(args, seed=42){
  nthreads = args$nthreads
  set.seed(seed)
  
  # read embeddings
  data <- read_csv(args$visualize.csv.gz)
  var_par <- as.data.frame(read_csv(args$simulate_var.csv.gz))
  # read SCE to get batch/celltype 
  sce <- read_sce(args$integrate.ad)
  batch <- sce$batch
  celltype <- sce$celltype
  rm(sce)
  
  val = sapply(unique(batch), function(b){
    h = var_par[var_par$batch==b,]
    l = apply(data[batch==b,],2, function(d) tapply(d, celltype[batch==b], var))
    l = l[complete.cases(l),]
    l = l[h$celltype,]
    h = h[complete.cases(h),3:ncol(h)]
    hv = rowSums(h)
    lv = rowSums(l)
    cor(hv, lv, method = "spearman")
              
  })
  return(mean(val))
}

variance_samplesize = function(args, seed=42){
  nthreads = args$nthreads
  set.seed(seed)
  
  # read embeddings
  data <- read_csv(args$visualize.csv.gz)
  # read SCE to get batch/celltype 
  sce <- read_sce(args$integrate.ad)
  batch <- sce$batch
  celltype <- sce$celltype
  rm(sce)
  
  val = sapply(unique(batch), function(b){
    h = tapply(which(batch==b), celltype[batch==b], length)
    l = rowSums(apply(data[batch==b,],2, function(d) tapply(d, celltype[batch==b], var)))
    h = h[complete.cases(h)]
    l = l[complete.cases(l)]
    1 - abs(cor(h, l, method = "spearman"))
  })
  return(mean(val))
}

library_size = function(args, seed=42){
  nthreads = args$nthreads
  B = args$B
  n = args$N
  set.seed(seed)
  
  # read embeddings
  data <- read_csv(args$visualize.csv.gz)
  # read SCE to get batch/celltype 
  sce <- read_sce(args$integrate.ad)
  batch <- sce$batch
  celltype <- sce$celltype
  library_size <- sce$nCount_RNA
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
        return(1 - tapply(ids, celltype[ids], function(idd) dcor(library_size[idd], data[idd,])$dcor))
      }else{
        return(NA)
      }
    }),na.rm=T)
    return(res)
  },mc.cores = nthreads)
  val = unlist(val)
  
  return(mean(val))
}

zero_proportion = function(args, seed=42){
  nthreads = args$nthreads
  B = args$B
  n = args$N
  set.seed(seed)
  
  # read embeddings
  data <- read_csv(args$visualize.csv.gz)
  # read SCE to get batch/celltype 
  sce <- read_sce(args$integrate.ad)
  batch <- sce$batch
  celltype <- sce$celltype
  mat <- counts(sce)
  zp <- colMeans(mat == 0)
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
        return(1 - tapply(ids, celltype[ids], function(idd) Rfast::dcor(zp[idd], data[idd,])$dcor))
      }else{
        return(NA)
      }
    }),na.rm=T)
    return(res)
  },mc.cores = nthreads)
  val = unlist(val)

  return(mean(val))
}
