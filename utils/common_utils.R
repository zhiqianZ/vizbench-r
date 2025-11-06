
suppressPackageStartupMessages({
  library(anndataR)
  library(rjson)
  library(readr)
})

read_sce <- function(f) read_h5ad(f, as = "SingleCellExperiment")
read_seurat <- function(f) read_h5ad(f, as = "Seurat")

write_ad <- function(x, file, verbose = TRUE) {
  if(verbose) message(paste("Converting", class(x), "-> AnnData."))
  if(typeof(x)!="environment"){
    x.ad <- as_AnnData(x)
  }else{
    x.ad <- x
  }
  if(verbose) message(paste0("Writing: ", file, "."))
  write_h5ad(x.ad, fn, mode = "w", compression = "gzip")
  if(verbose) message("Done.")
}

read_normmethod <- function(f) fromJSON(paste(readLines(f), 
                                              collapse=""))$normalize
