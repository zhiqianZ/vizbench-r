load_pkgs <- function() {
  library(scDesign3)
  library(SingleCellExperiment)
  library(Seurat)
  library(pbmcapply)
  library(parallel)
  library(readr)
}

.prepare_real_counts <- function(args) {
  sce <- read_sce(args$rawdata.ad)
  
  counts <- counts(sce)
  coldat <- colData(sce)
  ngenes <- args$ngenes
  
  meta <- data.frame(
    celltype = as.character(coldat$celltype),
    batch = as.character(coldat$batch)
  )
  
  ## filter genes detected in fewer than 10 cells
  filtered_genes <- which(Matrix::rowSums(counts != 0) < 10)
  if (length(filtered_genes) > 0) {
    counts <- counts[-filtered_genes, , drop = FALSE]
  }
  
  ## filter small batch-celltype combinations
  tab <- table(meta$celltype, meta$batch)
  small <- which(tab < 5 & tab > 0, arr.ind = TRUE)
  
  if (nrow(small) > 0) {
    small_pairs <- paste0(
      colnames(tab)[small[, "col"]],
      "_",
      rownames(tab)[small[, "row"]]
    )
    
    keep <- !paste0(meta$batch, "_", meta$celltype) %in% small_pairs
    
    counts <- counts[, keep, drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
    coldat <- coldat[keep, , drop = FALSE]
  }
  
  ## HVG selection
  seurat <- CreateSeuratObject(counts)
  seurat <- FindVariableFeatures(seurat, nfeatures = ngenes)
  hvg <- VariableFeatures(seurat)
  
  counts <- GetAssayData(seurat, layer = "counts")[hvg, , drop = FALSE]
  
  if (args$verbose) {
    message("Prepared real count matrix: ", paste0(dim(counts), collapse = ","))
  }
  
  list(
    counts = counts,
    coldat = coldat
  )
}

.fit_scdesign3 <- function(counts, coldat, args) {
  nthreads <- args$nthreads
  ncells <- args$ncells
  
  sce <- SingleCellExperiment(
    assay = list(counts = counts),
    colData = coldat
  )
  
  set.seed(123)
  
  if (args$verbose) {
    message("colData columns: ", paste0(colnames(colData(sce)), collapse = ", "))
  }
  
  if ("condition" %in% colnames(colData(sce))) {
    data <- construct_data(
      sce = sce,
      assay_use = "counts",
      celltype = "celltype",
      pseudotime = NULL,
      spatial = NULL,
      other_covariates = c("batch", "condition"),
      corr_by = "celltype",
      ncell = ncells
    )
  } else {
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
  }
  
  if (args$verbose) {
    message("Nonzero proportion in scDesign3 data: ", mean(data$count_mat > 0))
    message("scDesign3 count matrix dimension: ", paste0(dim(data$count_mat), collapse = ","))
  }
  
  message("Fitting marginal model")
  
  marginal <- fit_marginal(
    data = data,
    predictor = "gene",
    mu_formula = "celltype+batch",
    sigma_formula = "celltype",
    family_use = "nb",
    n_cores = nthreads,
    usebam = FALSE,
    parallelization = "mcmapply",
    trace = TRUE
  )
  
  message("Extracting parameters")
  
  para <- extract_para(
    sce = sce,
    marginal_list = marginal,
    n_cores = nthreads,
    family_use = "nb",
    new_covariate = data$newCovariate,
    data = data$dat,
    parallelization = "pbmcmapply"
  )
  
  list(
    sce = sce,
    data = data,
    marginal = marginal,
    para = para
  )
}

.make_parameter_tables <- function(seurat.obj, para) {
  batch <- seurat.obj$batch
  celltype <- seurat.obj$celltype
  uniq_batch <- unique(batch)
  
  idx_list <- setNames(
    lapply(uniq_batch, function(b) {
      by(which(batch == b), celltype[batch == b], \(id) id[1]) |> as.list()
    }),
    uniq_batch
  )
  
  mean_par <- lapply(names(idx_list), function(b) {
    idx <- idx_list[[b]]
    mat <- para$mean_mat[unlist(idx), , drop = FALSE]
    
    data.frame(
      batch = b,
      celltype = names(idx),
      mat,
      check.names = FALSE
    )
  })
  
  mean_par <- do.call(rbind, mean_par)
  mean_par <- mean_par[complete.cases(mean_par), ]
  
  var_par <- lapply(names(idx_list), function(b) {
    idx <- idx_list[[b]]
    
    mean_mat <- para$mean_mat[unlist(idx), , drop = FALSE]
    var_mat <- mean_mat^2 * para$sigma_mat[unlist(idx), , drop = FALSE]
    
    data.frame(
      batch = b,
      celltype = names(idx),
      var_mat,
      check.names = FALSE
    )
  })
  
  var_par <- do.call(rbind, var_par)
  var_par <- var_par[complete.cases(var_par), ]
  
  list(
    mean_par = mean_par,
    var_par = var_par
  )
}

.find_existing_parameter_files <- function(args) {
  mean_file_candidates <- c(
    args$simulate_mean.csv.gz,
    file.path(args$output_dir, paste0(args$name, "_simulate_mean.csv.gz"))
  )

  var_file_candidates <- c(
    args$simulate_var.csv.gz,
    file.path(args$output_dir, paste0(args$name, "_simulate_var.csv.gz"))
  )

  mean_file_candidates <- mean_file_candidates[
    !is.na(mean_file_candidates) & nzchar(mean_file_candidates)
  ]

  var_file_candidates <- var_file_candidates[
    !is.na(var_file_candidates) & nzchar(var_file_candidates)
  ]

  mean_file <- mean_file_candidates[file.exists(mean_file_candidates)][1]
  var_file <- var_file_candidates[file.exists(var_file_candidates)][1]

  if (is.na(mean_file) || is.na(var_file)) {
    return(NULL)
  }

  list(
    mean_file = mean_file,
    var_file = var_file
  )
}


scdesign3 <- function(args) {
  prepared <- .prepare_real_counts(args)
  
  counts <- prepared$counts
  coldat <- prepared$coldat
  
  fit <- .fit_scdesign3(
    counts = counts,
    coldat = coldat,
    args = args
  )
  
  message("Fitting copula")
  
  copula <- fit_copula(
    sce = fit$sce,
    assay_use = "counts",
    marginal_list = fit$marginal,
    family_use = "nb",
    copula = "gaussian",
    n_cores = args$nthreads,
    input_data = fit$data$dat,
    parallelization = "pbmcmapply"
  )
  
  rm(fit$marginal)
  gc()
  
  set.seed(123)
  
  message("Simulating new data")
  
  newcounts <- simu_new(
    sce = fit$sce,
    filtered_gene = fit$data$filtered_gene,
    mean_mat = fit$para$mean_mat,
    sigma_mat = fit$para$sigma_mat,
    zero_mat = fit$para$zero_mat,
    quantile_mat = NULL,
    copula_list = copula$copula_list,
    n_cores = args$nthreads,
    family_use = "nb",
    input_data = fit$data$dat,
    new_covariate = fit$data$newCovariate,
    important_feature = copula$important_feature,
    parallelization = "pbmcmapply"
  )
  
  seurat.obj <- CreateSeuratObject(
    newcounts,
    meta.data = fit$data$newCovariate
  )
  
  par_tabs <- .make_parameter_tables(
    seurat.obj = seurat.obj,
    para = fit$para
  )
  
  list(
    obj = seurat.obj,
    mean_par = par_tabs$mean_par,
    var_par = par_tabs$var_par
  )
}

real <- function(args) {
  prepared <- .prepare_real_counts(args)

  seurat.obj <- Seurat::CreateSeuratObject(
    prepared$counts,
    meta.data = prepared$coldat
  )

  par_files <- .find_existing_parameter_files(args)

  if (!is.null(par_files)) {
    message("Reading existing mean parameters: ", par_files$mean_file)
    message("Reading existing variance parameters: ", par_files$var_file)

    mean_par <- readr::read_csv(
      par_files$mean_file,
      show_col_types = FALSE
    )

    var_par <- readr::read_csv(
      par_files$var_file,
      show_col_types = FALSE
    )

  } else {
    message("Existing mean/variance parameter files not found.")
    message("Fitting scDesign3 marginal model to obtain parameters for real data.")

    fit <- .fit_scdesign3(
      counts = prepared$counts,
      coldat = prepared$coldat,
      args = args
    )

    par_tabs <- .make_parameter_tables(
      seurat.obj = seurat.obj,
      para = fit$para
    )

    mean_par <- par_tabs$mean_par
    var_par <- par_tabs$var_par
  }

  count_genes <- rownames(seurat.obj)
  mean_genes <- colnames(mean_par)[-(1:2)]
  var_genes <- colnames(var_par)[-(1:2)]

  if (!identical(count_genes, mean_genes)) {
    stop(
      "Gene mismatch between real counts and mean_par.\n",
      "This usually means --ngenes, filtering, or HVG selection differs between runs."
    )
  }

  if (!identical(count_genes, var_genes)) {
    stop(
      "Gene mismatch between real counts and var_par.\n",
      "This usually means --ngenes, filtering, or HVG selection differs between runs."
    )
  }

  list(
    obj = seurat.obj,
    mean_par = mean_par,
    var_par = var_par
  )
}






