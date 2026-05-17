load_pkgs <- function() {
  library(scDesign3)
  library(SingleCellExperiment)
  library(Seurat)
  library(pbmcapply)
  library(parallel)
  library(readr)
}

.align_parameter_genes <- function(par_df, count_genes, label_cols = c("batch", "celltype")) {
  par_genes <- setdiff(colnames(par_df), label_cols)

  if (!setequal(count_genes, par_genes)) {
    return(NULL)
  }

  par_df[, c(label_cols, count_genes), drop = FALSE]
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

.fit_marginal <- function(counts, coldat, args) {
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

.find_existing_parameter_files <- function(args, max_up = 5) {
  start_dir <- normalizePath(args$output_dir, mustWork = FALSE)

  if (!dir.exists(start_dir)) {
    start_dir <- dirname(start_dir)
  }

  ## Find the nearest ancestor named "simulate".
  ## Stop there. Do not continue beyond the dataset-level simulate directory.
  curr <- start_dir
  simulate_dir <- NULL

  for (i in seq_len(max_up)) {
    if (basename(curr) == "simulate") {
      simulate_dir <- curr
      break
    }

    parent <- dirname(curr)
    if (identical(parent, curr)) break
    curr <- parent
  }

  if (is.null(simulate_dir)) {
    message("No ancestor directory named 'simulate' found from: ", start_dir)
    return(NULL)
  }

  scdesign3_dir <- file.path(simulate_dir, "scdesign3")

  if (!dir.exists(scdesign3_dir)) {
    message("No sibling scdesign3 directory found at: ", scdesign3_dir)
    return(NULL)
  }

  message("Searching existing scDesign3 parameters under: ", scdesign3_dir)

  mean_files <- list.files(
    scdesign3_dir,
    pattern = "_simulate_mean\\.csv\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )

  var_files <- list.files(
    scdesign3_dir,
    pattern = "_simulate_var\\.csv\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )

  mean_files <- mean_files[file.exists(mean_files)]
  var_files <- var_files[file.exists(var_files)]

  if (length(mean_files) == 0 || length(var_files) == 0) {
    message("No mean/var files found under: ", scdesign3_dir)
    return(NULL)
  }

  mean_prefix <- sub("_simulate_mean\\.csv\\.gz$", "", basename(mean_files))
  var_prefix <- sub("_simulate_var\\.csv\\.gz$", "", basename(var_files))

  common_prefix <- intersect(mean_prefix, var_prefix)

  if (length(common_prefix) == 0) {
    message("Found mean/var files, but could not pair them by prefix.")
    return(NULL)
  }

  ## Prefer exact args$name match if available.
  exact_prefix <- common_prefix[common_prefix == args$name]

  if (length(exact_prefix) == 1) {
    chosen_prefix <- exact_prefix
  } else if (length(common_prefix) == 1) {
    chosen_prefix <- common_prefix
  } else {
    stop(
      "Multiple scDesign3 parameter file pairs found under:\n",
      scdesign3_dir,
      "\nCandidates:\n",
      paste(common_prefix, collapse = "\n"),
      "\nCannot safely choose one. Please ensure there is only one pair in this dataset's scdesign3 directory, or make the scdesign3 output prefix match args$name.",
      call. = FALSE
    )
  }

  mean_file <- mean_files[mean_prefix == chosen_prefix][1]
  var_file <- var_files[var_prefix == chosen_prefix][1]

  list(
    mean_file = normalizePath(mean_file),
    var_file = normalizePath(var_file)
  )
}


scdesign3 <- function(args) {
  prepared <- .prepare_real_counts(args)
  
  counts <- prepared$counts
  coldat <- prepared$coldat
  
  fit <- .fit_marginal(
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
  
  fit$marginal <- NULL
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

  count_genes <- rownames(seurat.obj)

  refit_parameters <- FALSE
  mean_par <- NULL
  var_par <- NULL

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

    mean_par_aligned <- .align_parameter_genes(mean_par, count_genes)
    var_par_aligned <- .align_parameter_genes(var_par, count_genes)

    if (is.null(mean_par_aligned) || is.null(var_par_aligned)) {
      message(
        "Gene set mismatch between real counts and existing mean/variance parameters. ",
        "Refitting scDesign3 marginal model."
      )
      refit_parameters <- TRUE
    } else {
      mean_par <- mean_par_aligned
      var_par <- var_par_aligned
    }

  } else {
    message("Existing mean/variance parameter files not found.")
    refit_parameters <- TRUE
  }

  if (refit_parameters) {
    message("Fitting scDesign3 marginal model to obtain parameters for real data.")

    fit <- .fit_marginal(
      counts = prepared$counts,
      coldat = prepared$coldat,
      args = args
    )

    par_tabs <- .make_parameter_tables(
      seurat.obj = seurat.obj,
      para = fit$para
    )

    mean_par <- .align_parameter_genes(par_tabs$mean_par, count_genes)
    var_par <- .align_parameter_genes(par_tabs$var_par, count_genes)

    if (is.null(mean_par)) {
      stop(
        "Gene set mismatch between real counts and refitted mean_par.\n",
        "Please check .prepare_real_counts(), .fit_marginal(), and .make_parameter_tables().",
        call. = FALSE
      )
    }

    if (is.null(var_par)) {
      stop(
        "Gene set mismatch between real counts and refitted var_par.\n",
        "Please check .prepare_real_counts(), .fit_marginal(), and .make_parameter_tables().",
        call. = FALSE
      )
    }
  }

  list(
    obj = seurat.obj,
    mean_par = mean_par,
    var_par = var_par
  )
}



