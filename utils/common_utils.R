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

Set_Threads_BLAS_OMP <- function(n = 255, init_python = TRUE, verbose = TRUE) {
  # Goal:
  # 1) Set *process env vars* that BOTH R and embedded Python (reticulate) can see.
  # 2) Set R-side thread pools (RhpcBLASctl, data.table).
  # 3) Optionally set Python-side pools (NumPy/BLAS via env; torch if present).
  #
  # IMPORTANT:
  # - For Python to truly honor env vars, set them BEFORE Python is initialized
  #   (i.e., before any reticulate::import() or py_config() that triggers init).
  #
  # Args:
  # - n: desired threads. If NULL, tries to infer a sane default from available cores.
  # - init_python: if TRUE and reticulate is installed, initialize Python and try to
  #   set torch threads if torch is available.
  # - verbose: print audits.

  # ---- infer threads (avoid hardcoding 255 unless you REALLY have it allocated) ----
  if (is.null(n)) {
    n <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
    if (is.na(n) || n <= 0) n <- 1L
  }
  n <- as.integer(n)

  # ---- helper to print env/thread status ----
  print_env <- function(header) {
    if (!verbose) return(invisible(NULL))
    message("\n== ", header, " ==")
    keys <- c(
      "OPENBLAS_NUM_THREADS","GOTO_NUM_THREADS","MKL_NUM_THREADS",
      "OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS",
      "OMP_DYNAMIC","MKL_DYNAMIC","OMP_PROC_BIND","KMP_AFFINITY",
      "KMP_BLOCKTIME","KMP_SETTINGS"
    )
    envs <- Sys.getenv(keys, unset = NA_character_)
    for (k in names(envs)) message(sprintf("  %s=%s", k, envs[[k]]))
  }

  print_r_threads <- function(header) {
    if (!verbose) return(invisible(NULL))
    message("\n== ", header, " ==")
    blas_path <- unname(extSoftVersion()[["BLAS"]])
    message(sprintf("BLAS: %s", if (is.na(blas_path)) "NA" else blas_path))

    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      bt <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NA_integer_)
      ot <- tryCatch(RhpcBLASctl::omp_get_max_threads(), error = function(e) NA_integer_)
      message(sprintf("R BLAS threads: %s", bt))
      message(sprintf("R OMP max threads: %s", ot))
    } else {
      message("RhpcBLASctl not installed (install.packages('RhpcBLASctl'))")
    }

    if (requireNamespace("data.table", quietly = TRUE)) {
      dtt <- tryCatch(data.table::getDTthreads(), error = function(e) NA_integer_)
      message(sprintf("data.table threads: %s", dtt))
    } else {
      message("data.table not installed")
    }
  }

  print_python_threads <- function(header) {
    if (!verbose) return(invisible(NULL))
    message("\n== ", header, " ==")
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      message("reticulate not installed; cannot audit Python")
      return(invisible(NULL))
    }
    # If Python hasn't been initialized, this will init it; only call if user asked.
    if (!reticulate::py_available(initialize = TRUE)) {
      message("Python not available in this session")
      return(invisible(NULL))
    }
    reticulate::py_run_string("
import os
print('Python sees env:')
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
          'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','OMP_DYNAMIC','MKL_DYNAMIC']:
    print(f'  {k}={os.getenv(k)}')
try:
    import numpy as np
    print('numpy:', np.__version__)
except Exception as e:
    print('numpy not available:', e)
try:
    from threadpoolctl import threadpool_info
    info = threadpool_info()
    print('threadpoolctl:', info)
except Exception as e:
    print('threadpoolctl not available:', e)
try:
    import torch
    print('torch:', torch.__version__)
    print('torch.get_num_threads():', torch.get_num_threads())
except Exception as e:
    print('torch not available:', e)
")
  }

  # ---- audit BEFORE ----
  print_env("Before: process env")
  print_r_threads("Before: R thread pools")

  # ---- 1) Set process env vars FIRST (shared by R + embedded Python) ----
  Sys.setenv(
    OPENBLAS_NUM_THREADS     = as.character(n),
    GOTO_NUM_THREADS         = as.character(n),
    MKL_NUM_THREADS          = as.character(n),
    OMP_NUM_THREADS          = as.character(n),
    VECLIB_MAXIMUM_THREADS   = as.character(n),  # macOS Accelerate
    NUMEXPR_NUM_THREADS      = as.character(n),
    OMP_DYNAMIC              = "FALSE",
    MKL_DYNAMIC              = "FALSE"
    # You can optionally add binding if you know what you're doing:
    # OMP_PROC_BIND          = "TRUE",
    # KMP_AFFINITY           = "granularity=fine,compact,1,0"
  )

  # ---- 2) Set R library thread pools (best-effort) ----
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    tryCatch(RhpcBLASctl::blas_set_num_threads(n), error = function(e) NULL)
    tryCatch(RhpcBLASctl::omp_set_num_threads(n),  error = function(e) NULL)
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    tryCatch(data.table::setDTthreads(n), error = function(e) NULL)
  }

  # ---- 3) Optional: initialize Python and set Python-side pools where possible ----
  if (init_python && requireNamespace("reticulate", quietly = TRUE)) {
    # Initializing Python AFTER setting env vars is key.
    if (reticulate::py_available(initialize = TRUE)) {
      # Some libs (e.g., torch) provide explicit setters.
      reticulate::py_run_string(sprintf("
try:
    import torch
    torch.set_num_threads(%d)
    torch.set_num_interop_threads(%d)
except Exception:
    pass
", n, max(1L, min(n, 8L))))  # interop threads often should be smaller
    }
  }

  # ---- audit AFTER ----
  print_env("After: process env")
  print_r_threads("After: R thread pools")
  if (init_python) print_python_threads("After: Python (reticulate) view")

  invisible(list(threads = n))
}
