#!/usr/bin/env Rscript

## run_vizbench.R
##
## Unified entry point for vizbench-r.
##
## Design:
##   --what    defines the pipeline stage, e.g. simulate, normalize,
##             integrate_normalize, visualize_count, metric_normalize.
##   --flavour defines the method/function to run. The function name must
##             exist in utils/<stage>_utils.R and must accept one argument:
##             args, a named list of parsed command-line arguments.
##
## Example:
##   Rscript run_vizbench.R \
##     --what simulate \
##     --flavour scdesign3 \
##     --rawdata.ad input/raw.h5ad \
##     --name dataset1 \
##     --output_dir output/
##
##   Rscript run_vizbench.R \
##     --what simulate \
##     --flavour real \
##     --rawdata.ad input/raw.h5ad \
##     --simulate_mean.csv.gz output/dataset1_simulate_mean.csv.gz \
##     --simulate_var.csv.gz output/dataset1_simulate_var.csv.gz \
##     --name dataset1_real \
##     --output_dir output/

suppressPackageStartupMessages({
  library(argparse)
  library(reticulate)
  library(jsonlite)
  library(readr)
})

## -----------------------------------------------------------------------------
## Helper functions
## -----------------------------------------------------------------------------

get_stage <- function(what) {
  sub("_.*$", "", what)
}

get_branch <- function(what) {
  if (!grepl("_", what)) {
    return(NA_character_)
  }
  sub("^.*_", "", what)
}

script_dir <- function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", cargs, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
  }

  ## Fallback for interactive/debug usage.
  return(getwd())
}

source_or_quit <- function(path, python = FALSE) {
  if (!file.exists(path)) {
    message("Helper code not found: ", path)
    quit("no", status = 1)
  }

  message("Sourcing: ", path)

  if (python) {
    reticulate::source_python(path, convert = FALSE)
  } else {
    source(path)
  }
}

require_arg <- function(args, name, stage_context = NULL) {
  value <- args[[name]]

  if (is.null(value) || is.na(value) || !nzchar(as.character(value))) {
    msg <- paste0("Missing required argument: --", name)
    if (!is.null(stage_context)) {
      msg <- paste0(msg, " for ", stage_context)
    }
    stop(msg, call. = FALSE)
  }

  invisible(value)
}

write_anndata_or_seurat <- function(x, path, verbose = TRUE) {
  if (verbose) {
    message("Writing AnnData/Seurat object: ", path)
  }

  if (typeof(x) != "environment") {
    write_seurat_ad(x, path)
  } else {
    x$write_h5ad(path, compression = "gzip")
  }
}

write_json <- function(x, path) {
  write(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE), path)
}

write_csv_gz <- function(x, path) {
  readr::write_csv(as.data.frame(x), file = path)
}

## -----------------------------------------------------------------------------
## Arguments
## -----------------------------------------------------------------------------

parser <- ArgumentParser(description = "vizbench-r benchmarking entry point")

parser$add_argument(
  "--what",
  choices = c(
    "rawdata",
    "simulate",
    "normalize",
    "integrate_normalize",
    "visualize_normalize",
    "integrate_count",
    "visualize_count",
    "metric_normalize",
    "metric_count"
  ),
  required = TRUE,
  help = paste(
    "Pipeline stage:",
    "rawdata, simulate, normalize, integrate_*, visualize_*, metric_*"
  )
)

parser$add_argument(
  "--flavour",
  choices = c(
    ## raw data
    "mouse_pancreas",
    "human_IFALD_liver",
    "human_atheroma",
    "human_covid_blood",
    "human_glaucoma_pbmc",
    "human_colorectal_liver",
    "human_prefrontal_cortex",
    "human_lung_atlas",
    "human_liver_atlas",
    "human_pbmc",
    "macaque_retina_fovea",
    "mouse_cortex",
    "mouse_intestine",
    "mouse_lung",
    "human_pancreas",
    "human_lung",
    "human_liver",

    ## simulate
    "scdesign3",
    "real",

    ## normalize
    "log1pCP10k",
    "log1pCPM",
    "sctransform",
    "log1pPF",
    "PFlog1pPF",
    "log1pCPMedian",
    "sanity",

    ## integrate
    "harmony",
    "fastMNN",
    "SeuratCCA",
    "SeuratRPCA",
    "LIGER",
    "scVI",
    "harmony-integrateRigor",
    "fastMNN-integrateRigor",
    "SeuratCCA-integrateRigor",
    "SeuratRPCA-integrateRigor",
    "LIGER-integrateRigor",
    "scVI-integrateRigor",

    ## visualize
    "SeuratUMAP",
    "scanpyUMAP",
    "BHtSNE",
    "FItSNE",
    "densMAP",
    "denSNE",
    "PHATE",
    "graphFA",
    "BHtSNE-scDEED",
    "SeuratUMAP-scDEED",
    "FItSNE-scDEED",
    "scanpyUMAP-scDEED",
    "densMAP-scDEED",

    ## metrics
    "celltype_shape",
    "batch_mixture",
    "batch_mixture_condition",
    "distance_preservation",
    "variance_preservation",
    "variance_samplesize",
    "library_size",
    "zero_proportion",
    "celltype_separation",
    "cLISI",
    "silhouette"
  ),
  required = TRUE,
  help = "Method/function to run. Must match a function name in utils/<stage>_utils.R."
)

## general benchmark parameters
parser$add_argument("--ncells", type = "integer", default = 1000000,
                    help = "Number of cells to simulate")
parser$add_argument("--ngenes", type = "integer", default = 3000,
                    help = "Number of genes to use/simulate")
parser$add_argument("--nthreads", type = "integer", default = 10,
                    help = "Number of threads")
parser$add_argument("--npcs", type = "integer", default = 50,
                    help = "Number of PCs")
parser$add_argument("--nhvgs", type = "integer", default = 2000,
                    help = "Number of HVGs")
parser$add_argument("--B", type = "integer", default = 100,
                    help = "Number of resampling replicates for metrics")
parser$add_argument("--N", type = "integer", default = 10000,
                    help = "Resampling sample size for metrics")

## use character instead of logical to avoid argparse TRUE/FALSE edge cases
parser$add_argument("--verbose", type = "character", default = "TRUE",
                    help = "TRUE/FALSE: whether to print progress messages")
parser$add_argument("--output_dir", "-o", dest = "output_dir",
                    type = "character", default = getwd(),
                    help = "Output directory where files will be saved")
parser$add_argument("--name", "-n", dest = "name",
                    type = "character", required = TRUE,
                    help = "Dataset/run name used as output file prefix")

## adaptive parameters tuning for scDEED and IntegrateRigor
parser$add_argument(
  "--parameters",
  type = "character",
  default = "{}",
  help = "JSON string specifying parameter ranges for scDEED and IntegrateRigor, e.g., --parameters '{\"theta\":[2,4,8],\"nclust\":[80,100,120]}'"
)

## stage input/output arguments
parser$add_argument("--rawdata.ad", type = "character",
                    help = "Input raw AnnData file")
parser$add_argument("--simulate.ad", type = "character",
                    help = "Input simulated/real-count AnnData file")
parser$add_argument("--simulate_mean.csv.gz", type = "character",
                    help = "CSV file containing batch-celltype gene-wise mean parameters")
parser$add_argument("--simulate_var.csv.gz", type = "character",
                    help = "CSV file containing batch-celltype gene-wise variance parameters")
parser$add_argument("--normalize.ad", type = "character",
                    help = "Input normalized AnnData file")
parser$add_argument("--normalize.json", type = "character",
                    help = "JSON file storing normalization method")
parser$add_argument("--integrate_normalize.ad", type = "character",
                    help = "Integrated AnnData from normalized branch")
parser$add_argument("--visualize_normalize.csv.gz", type = "character",
                    help = "2D embedding CSV from normalized branch")
parser$add_argument("--integrate_count.ad", type = "character",
                    help = "Integrated AnnData from count branch")
parser$add_argument("--visualize_count.csv.gz", type = "character",
                    help = "2D embedding CSV from count branch")

## Python / reticulate
parser$add_argument("--py_path", type = "character", default = "/usr/bin/python3",
                    help = "Python path for reticulate")
parser$add_argument("--scvi_conda", type = "character", default = "/usr/bin/python3",
                    help = "Optional scVI conda/Python environment identifier")

args <- parser$parse_args()

args$verbose <- tolower(args$verbose) %in% c("true", "t", "1", "yes", "y")

stage <- get_stage(args$what)
branch <- get_branch(args$what)
run_dir <- script_dir()

if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
}

## -----------------------------------------------------------------------------
## Logging
## -----------------------------------------------------------------------------

message("Full command: ", paste(commandArgs(), collapse = " "))
message("Selected stage: ", args$what)
message("Base stage: ", stage)
message("Branch: ", ifelse(is.na(branch), "none", branch))
message("Routine selected: ", args$flavour)
message("Additional parameters for tuning: ", args$parameters)
message("Name: ", args$name)
message("Output directory: ", args$output_dir)
message("Verbose: ", args$verbose)
message("Run directory: ", run_dir)
message("libPaths: ", paste(.libPaths(), collapse = ";"))

info <- Sys.info()
message("Sys.info: ", paste0(names(info), "=", info, collapse = ";"))

## -----------------------------------------------------------------------------
## Infer branch-specific inputs
## -----------------------------------------------------------------------------

if (stage == "visualize") {
  integrate_arg <- paste0("integrate_", branch, ".ad")
  args$integrate.ad <- args[[integrate_arg]]
  require_arg(args, "integrate.ad", args$what)
}

if (stage == "metric") {
  visualize_arg <- paste0("visualize_", branch, ".csv.gz")
  integrate_arg <- paste0("integrate_", branch, ".ad")

  args$visualize.csv.gz <- args[[visualize_arg]]
  args$integrate.ad <- args[[integrate_arg]]

  require_arg(args, "visualize.csv.gz", args$what)
  require_arg(args, "integrate.ad", args$what)
}

## Minimal input checks for common stages
if (args$what == "simulate") {
  require_arg(args, "rawdata.ad", "simulate")
}

if (args$what == "normalize") {
  require_arg(args, "simulate.ad", "normalize")
}

if (stage == "integrate") {
  if (branch == "normalize") {
    require_arg(args, "normalize.ad", args$what)
  }
  if (branch == "count") {
    require_arg(args, "simulate.ad", args$what)
  }
}

## -----------------------------------------------------------------------------
## Python setup
## -----------------------------------------------------------------------------

if (!is.null(args$py_path) && !is.na(args$py_path) && nzchar(args$py_path)) {
  reticulate::use_python(args$py_path, required = FALSE)
}

options(future.globals.maxSize = 10^20)

## -----------------------------------------------------------------------------
## Source helper files
## -----------------------------------------------------------------------------

common_helper <- file.path(run_dir, "utils", "common_utils.R")
source_or_quit(common_helper)

if (exists("Set_Threads_BLAS_OMP", mode = "function")) {
  Set_Threads_BLAS_OMP()
}

## source Python helper only for stages that use Python helper code
if (stage %in% c("normalize", "visualize")) {
  py_helper <- file.path(run_dir, "utils", paste0(stage, "_utils.py"))
  source_or_quit(py_helper, python = TRUE)
}

if (args$flavour == "FItSNE") {
  fitsne_path <- "/FIt-SNE/fast_tsne.R"
  if (file.exists(fitsne_path)) {
    source(fitsne_path, chdir = TRUE)
  } else {
    warning("FItSNE helper not found at: ", fitsne_path)
  }
}

stage_helper <- file.path(run_dir, "utils", paste0(stage, "_utils.R"))
source_or_quit(stage_helper)

## -----------------------------------------------------------------------------
## Load stage-specific packages
## -----------------------------------------------------------------------------

if (!exists("load_pkgs", mode = "function")) {
  stop("The helper file must define load_pkgs().", call. = FALSE)
}

suppressPackageStartupMessages(load_pkgs())

## -----------------------------------------------------------------------------
## Run selected method
## -----------------------------------------------------------------------------

fun <- tryCatch(
  get(args$flavour, mode = "function"),
  error = function(e) e
)

if (inherits(fun, "error")) {
  message("Unimplemented functionality: ", args$flavour)
  quit("no", status = 1)
}

message("Running: ", args$flavour)

x <- fun(as.list(args))

mean_par <- NULL
var_par <- NULL

if (args$what == "simulate") {
  if (!is.list(x) || !all(c("obj", "mean_par", "var_par") %in% names(x))) {
    stop(
      "simulate flavour must return a list with elements: obj, mean_par, var_par.",
      call. = FALSE
    )
  }

  mean_par <- x$mean_par
  var_par <- x$var_par
  x <- x$obj
  message("Object dimensions: ", paste(dim(x), collapse = " x "))
  if ("celltype" %in% colnames(x@meta.data)) {
    message(
    "Celltype summary:\n",
    paste(capture.output(table(x$celltype, useNA = "ifany")), collapse = "\n")
    )
  } else {
    message("Metadata column 'celltype' not found; skipping celltype table.")
  }
  if ("batch" %in% colnames(x@meta.data)) {
    message(
      "Batch summary:\n",
      paste(capture.output(table(x$batch, useNA = "ifany")), collapse = "\n")
    )
  } else {
    message("Metadata column 'batch' not found; skipping batch table.")
  }
  if ("condition" %in% colnames(x@meta.data)) {
    message(
      "Condition summary:\n",
      paste(capture.output(table(x$condition, useNA = "ifany")), collapse = "\n")
  )
  } else {
    message("Metadata column 'condition' not found; skipping condition table.")
  }
}

message("Done running: ", args$flavour)

## -----------------------------------------------------------------------------
## Save outputs
## -----------------------------------------------------------------------------

## AnnData-like outputs
if (stage %in% c("rawdata", "simulate", "normalize", "integrate")) {
  ad_path <- file.path(args$output_dir, paste0(args$name, "_", args$what, ".ad"))
  write_anndata_or_seurat(x, ad_path, verbose = args$verbose)
}

## Method metadata
if (args$what == "normalize") {
  json_path <- file.path(args$output_dir, paste0(args$name, "_", args$what, ".json"))
  write_json(list(normalize = args$flavour), json_path)
}

if (stage == "integrate") {
  json_path <- file.path(args$output_dir, paste0(args$name, "_", args$what, ".json"))
  write_json(list(integrate = args$flavour), json_path)
}

## Simulation parameter outputs
## Simulation parameter and metadata outputs
if (args$what == "simulate") {
  prefix <- file.path(args$output_dir, paste0(args$name, "_", args$what))

  mean_path <- paste0(prefix, "_mean.csv.gz")
  var_path <- paste0(prefix, "_var.csv.gz")
  celltype_path <- paste0(prefix, "_celltype.csv.gz")
  batch_path <- paste0(prefix, "_batch.csv.gz")
  condition_path <- paste0(prefix, "_condition.csv.gz")

  if (args$verbose) {
    message("Writing mean parameters: ", mean_path)
  }
  write_csv_gz(mean_par, mean_path)

  if (args$verbose) {
    message("Writing variance parameters: ", var_path)
  }
  write_csv_gz(var_par, var_path)

  write_meta_col <- function(obj, col, path) {
    if (col %in% colnames(obj@meta.data)) {
      out <- data.frame(
        cell = colnames(obj),
        value = as.character(obj[[col, drop = TRUE]]),
        stringsAsFactors = FALSE
      )
    } else {
      out <- data.frame(
        cell = colnames(obj),
        value = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    if (args$verbose) {
      message("Writing ", col, " metadata: ", path)
    }

    readr::write_csv(out, path)
  }
  write_meta_col(x, "celltype", celltype_path)
  write_meta_col(x, "batch", batch_path)
  write_meta_col(x, "condition", condition_path)
}

## Visualization outputs
if (stage == "visualize") {
  csv_path <- file.path(args$output_dir, paste0(args$name, "_", args$what, ".csv.gz"))

  if (args$verbose) {
    message("Writing visualization embedding: ", csv_path)
  }

  if (!args$flavour %in% c("scanpyUMAP", "graphFA")) {
    write_csv_gz(x, csv_path)
  } else {
    x_R <- reticulate::py_to_r(x)
    con <- gzfile(csv_path, open = "wt")
    on.exit(close(con), add = TRUE)
    write.csv(x_R, con, row.names = FALSE)
  }
}

## Metric outputs
if (stage == "metric") {
  json_path <- file.path(args$output_dir, paste0(args$name, "_", args$what, ".json"))

  if (args$verbose) {
    message("Writing metric: ", json_path)
  }

  write_json(list(value = x), json_path)
}

message("All outputs saved successfully.")
