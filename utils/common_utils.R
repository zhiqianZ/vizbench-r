suppressPackageStartupMessages({
  library(anndataR)
  library(rjson)
  library(readr)
})

read_sce <- function(f) read_h5ad(f, as = "SingleCellExperiment")
read_seurat <- function(f) read_h5ad(f, as = "Seurat")

write_seurat_ad <- function(x, file, verbose = TRUE) {
  if(verbose) message(paste("Converting", class(x), "-> AnnData."))
  x.ad <- as_AnnData(x)
  if(verbose) message(paste0("Writing: ", file, "."))
  write_h5ad(x.ad, fn, mode = "w", compression = "gzip")
  if(verbose) message("Done.")
}

read_normmethod <- function(f) fromJSON(paste(readLines(f), 
                                              collapse=""))$normalize

read_hvgs <- function(f) fromJSON(paste(readLines(f), 
                                              collapse=""))$hvgs

Set_Threads_BLAS_OMP <- function(){
  ## == OpenMP / threading env ==
  envs <- Sys.getenv(c(
    "OPENBLAS_NUM_THREADS","GOTO_NUM_THREADS","MKL_NUM_THREADS",
    "OMP_NUM_THREADS","OMP_DYNAMIC","MKL_DYNAMIC",
    "OMP_PROC_BIND","KMP_AFFINITY"
  ), unset = NA_character_)
  message("== OpenMP / threading env ==")
  for (k in names(envs)) message(sprintf("  %s=%s", k, envs[[k]]))
  
  ## == R: BLAS & thread counts ==
  message("\n== R: BLAS & thread counts ==")
  blas_path <- unname(extSoftVersion()[["BLAS"]])
  message(sprintf("BLAS: %s", if (is.na(blas_path)) "NA" else blas_path))
  
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    bt <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error=function(e) NA_integer_)
    ot <- tryCatch(RhpcBLASctl::omp_get_max_threads(), error=function(e) NA_integer_)
    message(sprintf("BLAS threads: %s", bt))
    message(sprintf("OMP max threads: %s", ot))
  } else {
    message("RhpcBLASctl not installed (install.packages('RhpcBLASctl'))")
  }
  
  if (requireNamespace("data.table", quietly = TRUE)) {
    dtt <- tryCatch(data.table::getDTthreads(), error=function(e) NA_integer_)
    message(sprintf("data.table threads: %s", dtt))
  }
  
  ## env vars first (so child libs pick them up)
  Sys.setenv(
    OPENBLAS_NUM_THREADS = "255",
    GOTO_NUM_THREADS     = "255",
    MKL_NUM_THREADS      = "255",
    OMP_NUM_THREADS      = "255",
    MKL_DYNAMIC          = "FALSE",
    OMP_DYNAMIC          = "FALSE"
  )
  
  ## set library thread pools
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(255)
    RhpcBLASctl::omp_set_num_threads(255)
    message(paste("BLAS threads:", RhpcBLASctl::blas_get_num_procs()))
    message(paste("OMP max threads:", RhpcBLASctl::omp_get_max_threads()))
  } else {
    message("RhpcBLASctl not installed; cannot set BLAS/OMP pools from R")
  }
  
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::setDTthreads(255)
    message(paste("data.table threads:", data.table::getDTthreads()))
  }
  
  ## show final env (audit)
  for (k in c("OPENBLAS_NUM_THREADS","GOTO_NUM_THREADS","MKL_NUM_THREADS",
              "OMP_NUM_THREADS","MKL_DYNAMIC","OMP_DYNAMIC")) {
    message(paste0(k, "=", Sys.getenv(k, unset = "NA")))
  }
  message("Updated")
  ## == OpenMP / threading env ==
  envs <- Sys.getenv(c(
    "OPENBLAS_NUM_THREADS","GOTO_NUM_THREADS","MKL_NUM_THREADS",
    "OMP_NUM_THREADS","OMP_DYNAMIC","MKL_DYNAMIC",
    "OMP_PROC_BIND","KMP_AFFINITY"
  ), unset = NA_character_)
  message("== OpenMP / threading env ==")
  for (k in names(envs)) message(sprintf("  %s=%s", k, envs[[k]]))
    
  ## == R: BLAS & thread counts ==
  message("\n== R: BLAS & thread counts ==")
  blas_path <- unname(extSoftVersion()[["BLAS"]])
  message(sprintf("BLAS: %s", if (is.na(blas_path)) "NA" else blas_path))
  
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    bt <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error=function(e) NA_integer_)
    ot <- tryCatch(RhpcBLASctl::omp_get_max_threads(), error=function(e) NA_integer_)
    message(sprintf("BLAS threads: %s", bt))
    message(sprintf("OMP max threads: %s", ot))
  } else {
    message("RhpcBLASctl not installed (install.packages('RhpcBLASctl'))")
  }
  
  if (requireNamespace("data.table", quietly = TRUE)) {
    dtt <- tryCatch(data.table::getDTthreads(), error=function(e) NA_integer_)
    message(sprintf("data.table threads: %s", dtt))
  }
  
}
