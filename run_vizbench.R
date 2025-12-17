#!/usr/bin/env Rscript

## Originally sketched by Izaskun Mallona
## modified/populated by Mark Robinson
## various chunks of code from Zhiqian Zhai and Qingyang Wang will be linked
## from original repo: https://github.com/zhiqianZ/Benchmark-Normalization-Integration-Visualization

## Usage:
## to do a system call to showcase how to call other subscripts
##    Rscript thisfilename.R --what type --flavour specific-module

library(argparse)
library(reticulate)

parser <- ArgumentParser(description = "Benchmarking entrypoint")

# define arguments
parser$add_argument("--what", 
                    choices = c("rawdata", "simulate", "normalize", 
                                "integratenorm", "integrateraw", "visualize", "metric"),
                    required = TRUE, 
                    help = "Module type: rawdata, simulate, normalize, integrate, vizualize, metric")

#s <- switch(args$what, rawdata = c("mouse_pancreas"),
#                       simulation = c("scdesign3"))
# TODO: add subparser?

parser$add_argument("--flavour", 
                    choices = c("mouse_pancreas",                                                                                    # raw data
                                "scdesign3",                                                                                         # simulate
                                "log1pCP10k", "log1pCPM", "sctransform", "log1pPF", "PFlog1pPF", "log1pCPMedian", "sanity",          # normalize
                                "harmony", "fastMNN", "SeuratCCA", "SeuratRPCA", "LIGER", "scVI",                                    # integrate
                                "SeuratUMAP", "scanpyUMAP", "BHtSNE", "FItSNE", "densMAP", "denSNE", "PHATE", "graphFA",             # visualize
                                "celltype_shape", "batch_mixture", "distance_preservation", "variance_preservation",
                                "variance_samplesize", "library_size", "zero_proportion", "celltype_separation"),                    # metrics
                    required = TRUE, 
                    help = "Module to run: name depends on the 'what'")

parser$add_argument("--ncells", type = "integer", default = 1000000,
                    help = "number of cells to simulate")

parser$add_argument("--ngenes", type = "integer", default = 2000,
                    help = "number of genes to simulate")

parser$add_argument("--nthreads", type = "integer", default = 10,
                    help = "numer of threads used for simulating and benchmarking methods")

parser$add_argument("--npcs", type = "integer", default = 50,
                    help = "number of pcs used")

parser$add_argument("--nhvgs", type = "integer", default = 2000,
                    help = "number of hvgs used")

parser$add_argument("--B", type = "integer", default = 100,
                    help = "Resampling replicate times for calcualting metrics")

parser$add_argument("--N", type = "integer", default = 10000,
                    help = "Resampling sample size for calcualting metrics")

parser$add_argument("--use_simulation", type = "logical", default = TRUE,
                    help = "whether used the simulated datasets to benchmark")

parser$add_argument("--verbose", type = "logical", default = TRUE,
                    help = "TRUE/FALSE as to whether to write progress to stdout")

parser$add_argument("--output_dir", "-o", dest="output_dir", type="character",
                    help="output directory where files will be saved", default=getwd(),
                    required = TRUE)

parser$add_argument("--name", "-n", dest="name", type="character", required = TRUE,
                    help="name of the dataset")

parser$add_argument('--rawdata.ad',
                    type="character",
                    help='gz-compressed H5 file containing (raw) data as AnnData')

parser$add_argument('--simulate.ad',
                    type="character",
                    help='gz-compressed H5 file containing (simulated) data as AnnData')

parser$add_argument('--simulate_mean.csv.gz',
                    type="character",
                    help='gz-compressed CSV file containing per batch and cell-type gene-wise true mean parameters')

parser$add_argument('--simulate_var.csv.gz',
                    type="character",
                    help='gz-compressed CSV file containing per batch and cell-type gene-wise true variance parameters')

parser$add_argument('--normalize.ad',
                    type="character",
                    help='gz-compressed H5 file containing (normalized) data as AnnData')

parser$add_argument('--normalize.json',
                    type="character",
                    help='JSON file containing name of normalization method')

parser$add_argument('--integrate.ad',
                    type="character",
                    help='gz-compressed H5 file containing (integrated_from_nrom) data as AnnData')

parser$add_argument('--visualize.csv.gz',
                    type="character",
                    help='gz-compressed CSV file containing embeddings')

parser$add_argument('--py_path', 
                    type="character",
                    default="/usr/bin/python3",
                    help='the path of the Python for the reticulate package to use')

parser$add_argument('--scvi_conda', 
                    type="character",
                    default="/usr/bin/python3",
                    help='the path of the Python for the reticulate package to use')

# send details to be logged
args <- parser$parse_args()
message("Full command: ", paste0(commandArgs(), collapse = " "))
message("Selected category: ", args$what)
message("Routine selected: ", args$flavour)
message("Additional parameters: ", args$params)
message("name: ", args$name)
message("Verbose: ", args$verbose)

# infer the current directory (useful for debugging)
cargs <- commandArgs(trailingOnly = FALSE)
m <- grep("--file=", cargs)
run_dir <- dirname( gsub("--file=","",cargs[[m]]) )
message("location: ", run_dir)
message("libPaths: ", paste0(.libPaths(),collapse=";"))
info <- Sys.info()
message("info: ", paste0(names(info),"=",info,collapse=";"))


if(args$what %in% c("integratenorm", "integrateraw")){
  args$what = "integrate"
}

options(future.globals.maxSize= 10^20)
Set_Threads_BLAS_OMP()

# source common helper functions
helpers <- file.path(run_dir, "utils", "common_utils.R")
if( file.exists(helpers) ) {
  message("Sourcing .. ", helpers)
  source(helpers)
} else {
  message(paste0("Helper code in ", helpers, " not found. Exiting."))
  quit("no", status = 1)
}

use_python(args$py_path)

# source normalization python helper functions
if (args$what %in% c("normalize", "visualize")) {
  helpers <- file.path(run_dir, "utils", paste0(args$what, "_utils.py"))
  if( file.exists(helpers) ) {
    message("Sourcing .. ", helpers)
    source_python(helpers, convert = FALSE)
  } else {
    message(paste0("Helper code in ", helpers, " not found. Exiting."))
    quit("no", status = 1)
  }
}
if(args$flavour == "FItSNE"){
  source("/FIt-SNE/fast_tsne.R",chdir=T)
}

# source stage-specific helper functions (n.b.: according to args$what)
helpers <- file.path(run_dir, "utils", paste0(args$what, "_utils.R"))
if( file.exists(helpers) ) {
    message("Sourcing .. ", helpers)
    source(helpers)
} else {
    message(paste0("Helper code in ", helpers, " not found. Exiting."))
    quit("no", status = 1)
} 

# load packages
suppressPackageStartupMessages(load_pkgs())

# check if implemented: throw error if not; run if so
# n.b.: args$flavour defines what 'main' function to call
fun <- tryCatch(obj <- get(args$flavour), error = function(e) e)
if ( !("error" %in% class(fun)) ) {
    x <- fun(as.list(args)) # execute function 
    if(args$what == "simulate"){
      mean_par <- x$mean_par
      var_par <- x$var_par
      x <- x$obj
    }
  message("done running")
} else {
    message('Unimplemented functionality. Exiting.\n') # throw error?
    quit("no", status = 1)
}

# write to AnnData via anndataR
if (args$what %in% c("rawdata", "simulate", "normalize", "integrate")) {
  # here, always writing data files as AD (HDF5)
  fn <- file.path(args$output_dir, paste0(args$name,"_",args$what, ".ad"))
  if(typeof(x)!="environment"){
    write_seurat_ad(x, fn)
  }else{
    if(args$verbose) message(paste0("Writing: ", fn, "."))
    x$write_h5ad(fn, compression = "gzip")
  }
} 
# write memento about normalization method
if(args$what == "normalize") {
  fn <- file.path(args$output_dir, paste0(args$name,"_", args$what, ".json"))
  write(toJSON(list(normalize=args$flavour)), fn)
}
                
if(args$what == "integrate"){
  fn <- file.path(args$output_dir, paste0(args$name, args$what, ".json"))
  write(toJSON(list(intgrate=args$flavour)), fn)
}
                
if(args$what == "simulate"){
  fn <- file.path(args$output_dir, paste0(args$name,"_",args$what, "_mean.csv.gz"))
  write_csv(as.data.frame(mean_par), file = fn)
  fn <- file.path(args$output_dir, paste0(args$name,"_",args$what, "_var.csv.gz"))
  write_csv(as.data.frame(var_par), file = fn)
}
                
if(args$what == "visualize") {
  # here, write embeddings to gzipped CSV file
  fn <- gzfile(file.path(args$output_dir, 
                         paste0(args$name, "_visualize", ".csv.gz")))
  if(!args$flavour %in% c("scanpyUMAP","graphFA")){
    write_csv(as.data.frame(x), file = fn)
  }else{
    x_R <- py_to_r(x)
    write.csv(x_R, fn, row.names = FALSE)
  }
} else if (args$what == "metric") {
  # 'x' is something here
  fn <- file.path(args$output_dir, paste0(args$name,"_",args$what, ".json"))
  write(toJSON(list(value=x)), fn)
}
