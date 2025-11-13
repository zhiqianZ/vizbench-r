
load_pkgs <- function() {
  library(Seurat)
}

log1pCP10k = function(args){
  message("Running log1pCP10k")
  seurat.obj <- read_seurat(args$simulate.ad)
  NormalizeData(seurat.obj, scale.factor = 10^4)
}

log1pCPM = function(args){
  message("Running log1pCPM")
  seurat.obj <- read_seurat(args$simulate.ad)
  NormalizeData(seurat.obj, scale.factor= 10^6)
}

sctransform = function(args){
  message("Running SCTransform")
  seurat.obj <- read_seurat(args$simulate.ad)
  nhvgs = args$nhvgs
  seurat.obj[["RNA"]] <- split(seurat.obj[["RNA"]], f = seurat.obj$batch)
  seurat.obj = SCTransform(seurat.obj, vst.flavor = "v2", verbose = FALSE, 
              variable.features.n = nhvgs, method = "glmGamPoi")
  
  seurat.obj = seurat.obj[rownames(seurat.obj[['SCT']]$scale.data),]
  return(seurat.obj)
}


Sanity = function(args){
  sanity_exec <- args$sanity_path
  input_mat <- args$input_mat
  input_genes <- args$input_genes 
  input_cells <- args$input_cells
  output_dir <- args$sanity_output
  num_threads = args$nthreads
  system2(file.path(sanity_exec, "bin/Sanity"), args = c(paste0("--file ", input_mat),
                                                                  paste0("--mtx_gene_name_file ", input_genes),
                                                                  paste0("--mtx_cell_name_file ", input_cells),
                                                                 paste0("--destination ", output_dir),
                                                                 paste0("--n_threads ", num_threads)))
  return("Done")
}









