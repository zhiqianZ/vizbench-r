load_pkgs <- function() {
  library(scDesign3)
  library(SingleCellExperiment)
  library(Seurat)
  library(pbmcapply)
  library(parallel)
}

scdesign3 <- function(args) {
  
  sce <- read_sce(args$rawdata.ad)
  counts = counts(sce)
  coldat = colData(sce)

  ncells = args$ncells
  ngenes = args$ngenes
  nthreads = args$nthreads
  use_simulation = args$use_simulation
	
  # QC
  batch = unique(coldat$batch)
  meta = data.frame(cbind(celltype=as.character(coldat$celltype),
                          batch=as.character(coldat$batch)))
  batch_filtered = colnames(table(meta))[which(table(meta)<5&table(meta)>0, arr.ind = T)[,2]]
  celltype_filtered = rownames(table(meta))[which(table(meta)<5&table(meta)>0, arr.ind = T)[,1]]
  filtered_genes = which(apply(counts,1,function(x) length(which(x!=0)))<10)
  if(length(filtered_genes)>0){
    counts = counts[-filtered_genes,]  
  }
  if(length(batch_filtered) > 0){
    filter = paste0(batch_filtered,"_", celltype_filtered)
    idx = which(paste0(meta$batch,"_", meta$celltype) %in% filter)
    counts = counts[, -idx]
    meta = meta[-idx, ]
    coldat = coldat[-idx,]
  }

  if(use_simulation == FALSE){ncells = ncol(counts)}
							   
  # HVG 
  seurat = CreateSeuratObject(counts)
  seurat = FindVariableFeatures(seurat,nfeatures = ngenes) 
  hvg = VariableFeatures(seurat)
  counts = counts[hvg,]
							   
  if(args$verbose) message(paste0(dim(counts), collapse=","))
  ## scDesign3
  sce <- SingleCellExperiment(assay = list(counts = counts), 
                              colData = coldat)

  set.seed(123)
  data <- construct_data(
    sce = sce,
    assay_use = "counts",
    celltype = "celltype",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = "batch",
    corr_by = "celltype",
    ncell = ncells
  )
  
  if(args$verbose) message(mean(data$count_mat>0))
  if(args$verbose) message(paste0(dim(data$count_mat), collapse=","))

  message("Fitting marginal")
  marginal <- fit_marginal(
    data = data,
    predictor = "gene",
    #mu_formula = "celltype+s(batch,bs='re')",
    mu_formula = "celltype+batch",
    sigma_formula = "celltype",
    family_use = "nb",
    n_cores = nthreads,
    usebam = FALSE,
    #parallelization = "mcmapply",
    parallelization = "pbmcmapply",
    trace = TRUE
  )
  para <- extract_para(
    sce = sce,
    marginal_list = marginal,
    n_cores = nthreads,
    family_use = "nb",
    new_covariate = data$newCovariate,
    #parallelization = "mcmapply",
    data = data$dat,
    parallelization = "pbmcmapply"
  )


  if(args$verbose) message(paste0(names(data), collapse=","))
  if(use_simulation == TRUE){
	  message("Fitting copula")
	  copula <- fit_copula(
		  sce = sce,
    	  assay_use = "counts",
    	  marginal_list = marginal,
    	  family_use = "nb",
    	  copula = "gaussian",
    	  n_cores = nthreads,
    	  input_data = data$dat,
		  parallelization = "pbmcmapply"
	  )
	  set.seed(123)
	  message("Simulating new data")						   
	  newcounts = simu_new(
    	  sce = sce,
    	  filtered_gene = data$filtered_gene,
    	  mean_mat = para$mean_mat,
    	  sigma_mat = para$sigma_mat,
    	  zero_mat = para$zero_mat,
    	  quantile_mat = NULL,
    	  copula_list = copula$copula_list,
    	  n_cores = nthreads,
    	  family_use = "nb",
    	  input_data = data$dat,
    	  new_covariate = data$newCovariate,
    	  important_feature = copula$important_feature,
    	  parallelization = "pbmcmapply"
  	  )
	  seurat.obj = CreateSeuratObject(newcounts, 
				  meta.data=data$newCovariate)
  }else{
	  seurat.obj = CreateSeuratObject(counts, 
				  meta.data=coldat)
  }
  batch <- coldat$batch
  celltype <- coldat$celltype
  uniq_batch <- unique(batch)

  #Build index list: one representative cell index per celltype in each batch
  idx_list <- setNames(
    lapply(uniq_batch, function(b) {
      by(which(batch == b), celltype[batch == b], \(id) id[1]) |> as.list()
    }),
    uniq_batch
  )

  # Compute mean and variance parameters
  parameters <- list()

  parameters$mean <- lapply(idx_list, function(idx) {
    mat <- para$mean_mat[unlist(idx), , drop = FALSE]
    rownames(mat) <- names(idx)
    mat
  })

  parameters$var <- lapply(idx_list, function(idx) {
    mean_mat <- para$mean_mat[unlist(idx), , drop = FALSE]
    var_mat <- mean_mat^2 * para$sigma_mat[unlist(idx), , drop = FALSE]
    rownames(var_mat) <- names(idx)
    var_mat
  })

						
  return(list(obj=seurat.obj, parameters=parameters))
}

