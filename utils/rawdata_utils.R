
load_pkgs <- function() {
  library(readr)
  library(SingleCellExperiment)
  library(Seurat)
  library(dplyr)
  library(GEOquery)
  library(R.utils)
  library(here)
  library(zellkonverter)
}

mouse_pancreas <- function(args) {
  ### Mouse pancreas GSM2230761 GSM2230762 (pilot)
  file = c("GSM2230761", "GSM2230762")
  data = lapply(1:2, function(x) {
    tmpdir <- tempdir()
    getGEOSuppFiles(file[x], baseDir = tmpdir)
    curr_data = read.csv(file.path(tmpdir, file[x],
				  paste0(file[x],"_mouse",x,"_umifm_counts.csv.gz")))
    return(curr_data)
  })
  
  label = lapply(data, function(obj) obj[,3])
  label = unlist(label)
  batch = rep(paste0("inDrop",1:2), unlist(lapply(data,nrow)))
  data = lapply(data,function(obj) {rownames(obj)=obj[,1]; obj[,-(1:3)]})
  count_mat = t(do.call("rbind", data))
  coldata = data.frame(celltype = label, batch = batch)
  rownames(coldata) = colnames(count_mat)
  sce_save = SingleCellExperiment(list(counts = count_mat), 
                                  colData = coldata)
  return(sce_save)
}

human_IFALD_liver <- function(args) {
  options(timeout = 7200)
  ### https://cellxgene.cziscience.com/collections/ff69f0ee-fef6-4895-9f48-6c64a68c8289 
  url = "https://datasets.cellxgene.cziscience.com/feca90bb-00df-4623-8398-1e3e6a90971d.h5ad"
  temp_h5ad <- tempfile(fileext = ".h5ad")
  download.file(url, destfile = temp_h5ad, mode = "wb")
  sce <- readH5AD(
    file = temp_h5ad, 
    X_name = "logcounts", 
    raw = TRUE
  )
  sce_raw <- as(SingleCellExperiment::altExp(sce, "raw"), "SingleCellExperiment")
  assayNames(sce_raw)[assayNames(sce_raw) == "X"] <- "counts"
  raw_values <- counts(sce_raw)
  is_integer <- all(raw_values == round(raw_values))
  
  if (is_integer) {
    message("The 'counts' assay contains only integers.")
  } else {
    warning("Warning: The 'counts' assay contains decimals. This data is likely already normalized.")
  }
  common_cells <- intersect(colnames(sce), colnames(sce_raw))
  sce_raw <- sce_raw[,common_cells]
  colData(sce_raw) <- colData(sce)[colnames(sce_raw),]
  colData(sce_raw)$celltype = colData(sce_raw)$cell_type
  colData(sce_raw)$batch = colData(sce_raw)$donor_id
  colData(sce_raw)$condition = colData(sce_raw)$disease
  cond <- as.character(colData(sce_raw)$condition)
  cond[cond != "normal"] <- "intestinal failure-associated liver disease"
  colData(sce_raw)$condition <- cond
  sce_raw = sce_raw[,sce_raw$celltype != "unknown"]
  message("celltype:\n", paste(capture.output(table(colData(sce_raw)$celltype)), collapse = "\n"))
  message("batch:\n", paste(capture.output(table(colData(sce_raw)$batch)), collapse = "\n"))
  message("condtion:\n", paste(capture.output(table(colData(sce_raw)$condition)), collapse = "\n"))
  file.remove(temp_h5ad)
  return(sce_raw)
}

human_atheroma <- function(args) {
  options(timeout = 7200)
  ### https://cellxgene.cziscience.com/collections/db70986c-7d91-49fe-a399-a4730be394ac
  ### therosclerotic plaque scRNAseq datasets
  url = "https://datasets.cellxgene.cziscience.com/3bab1c0b-d3e3-4a01-840f-d49a8284d989.h5ad"

  temp_h5ad <- tempfile(fileext = ".h5ad")
  download.file(url, destfile = temp_h5ad, mode = "wb")
  sce <- readH5AD(
    file = temp_h5ad, 
    X_name = "logcounts", 
    raw = TRUE
  )
  sce_raw <- as(SingleCellExperiment::altExp(sce, "raw"), "SingleCellExperiment")
  assayNames(sce_raw)[assayNames(sce_raw) == "X"] <- "counts"
  raw_values <- counts(sce_raw)
  is_integer <- all(raw_values == round(raw_values))
  
  if (is_integer) {
    message("The 'counts' assay contains only integers.")
  } else {
    warning("Warning: The 'counts' assay contains decimals. This data is likely already normalized.")
  }
  common_cells <- intersect(colnames(sce), colnames(sce_raw))
  sce_raw <- sce_raw[,common_cells]
  colData(sce_raw) <- colData(sce)[colnames(sce_raw),]
  colData(sce_raw)$celltype = colData(sce_raw)$cell_type
  colData(sce_raw)$batch = colData(sce_raw)$donor_id
  # colData(sce_raw)$condition = colData(sce_raw)$disease
  message("celltype:\n", paste(capture.output(table(colData(sce_raw)$celltype)), collapse = "\n"))
  message("batch:\n", paste(capture.output(table(colData(sce_raw)$batch)), collapse = "\n"))
  file.remove(temp_h5ad)
  return(sce_raw)
}


human_glaucoma_pbmc <- function(args) {
  options(timeout = 7200)
  ### https://cellxgene.cziscience.com/collections/de2cde16-c8d3-4a6d-80be-1be9e879aaca
  url = "https://datasets.cellxgene.cziscience.com/06932380-c04d-4d11-b4af-22de62b031b4.h5ad"

  temp_h5ad <- tempfile(fileext = ".h5ad")
  download.file(url, destfile = temp_h5ad, mode = "wb")
  sce <- readH5AD(
    file = temp_h5ad, 
    X_name = "logcounts", 
    raw = TRUE
  )
  sce_raw <- as(SingleCellExperiment::altExp(sce, "raw"), "SingleCellExperiment")
  assayNames(sce_raw)[assayNames(sce_raw) == "X"] <- "counts"
  raw_values <- counts(sce_raw)
  is_integer <- all(raw_values == round(raw_values))
  
  if (is_integer) {
    message("The 'counts' assay contains only integers.")
  } else {
    warning("Warning: The 'counts' assay contains decimals. This data is likely already normalized.")
  }
  common_cells <- intersect(colnames(sce), colnames(sce_raw))
  sce_raw <- sce_raw[,common_cells]
  colData(sce_raw) <- colData(sce)[colnames(sce_raw),]
  colData(sce_raw)$celltype = colData(sce_raw)$cell_type
  colData(sce_raw)$batch = colData(sce_raw)$donor_id
  colData(sce_raw)$condition = colData(sce_raw)$disease
  message("celltype:\n", paste(capture.output(table(colData(sce_raw)$celltype)), collapse = "\n"))
  message("batch:\n", paste(capture.output(table(colData(sce_raw)$batch)), collapse = "\n"))
  message("condtion:\n", paste(capture.output(table(colData(sce_raw)$condition)), collapse = "\n"))
  file.remove(temp_h5ad)
  return(sce_raw)
}
				 
human_covid_blood <- function(args) {
  options(timeout = 7200)
  ### https://cellxgene.cziscience.com/collections/8f126edf-5405-4731-8374-b5ce11f53e82
  url = "https://datasets.cellxgene.cziscience.com/687c09ff-731a-4e3d-ac07-4c29c33a6338.h5ad"
  temp_h5ad <- tempfile(fileext = ".h5ad")
  download.file(url, destfile = temp_h5ad, mode = "wb")
  sce <- readH5AD(
    file = temp_h5ad, 
    X_name = "logcounts", 
    raw = TRUE
  )
  sce_raw <- as(SingleCellExperiment::altExp(sce, "raw"), "SingleCellExperiment")
  assayNames(sce_raw)[assayNames(sce_raw) == "X"] <- "counts"
  raw_values <- counts(sce_raw)
  is_integer <- all(raw_values == round(raw_values))
  
  if (is_integer) {
    message("The 'counts' assay contains only integers.")
  } else {
    warning("Warning: The 'counts' assay contains decimals. This data is likely already normalized.")
  }
  common_cells <- intersect(colnames(sce), colnames(sce_raw))
  sce_raw <- sce_raw[,common_cells]
  colData(sce_raw) <- colData(sce)[colnames(sce_raw),]
  colData(sce_raw)$celltype = colData(sce_raw)$cell_type
  colData(sce_raw)$batch = colData(sce_raw)$donor_id
  colData(sce_raw)$condition = colData(sce_raw)$disease
  message("celltype:\n", paste(capture.output(table(colData(sce_raw)$celltype)), collapse = "\n"))
  message("batch:\n", paste(capture.output(table(colData(sce_raw)$batch)), collapse = "\n"))
  message("condtion:\n", paste(capture.output(table(colData(sce_raw)$condition)), collapse = "\n"))
  file.remove(temp_h5ad)
  return(sce_raw)
}


human_colorectal_liver <- function(args) {
  options(timeout = 7200)
  ### https://cellxgene.cziscience.com/collections/be679cb1-35f0-46c9-9a2d-30691862a54a
  url = "https://datasets.cellxgene.cziscience.com/b425976f-9d73-4388-95dd-e7cd0f8caca0.h5ad"
  temp_h5ad <- tempfile(fileext = ".h5ad")
  download.file(url, destfile = temp_h5ad, mode = "wb")
  sce <- readH5AD(
    file = temp_h5ad, 
    X_name = "logcounts", 
    raw = TRUE
  )
  sce_raw <- as(SingleCellExperiment::altExp(sce, "raw"), "SingleCellExperiment")
  assayNames(sce_raw)[assayNames(sce_raw) == "X"] <- "counts"
  raw_values <- counts(sce_raw)
  is_integer <- all(raw_values == round(raw_values))
  
  if (is_integer) {
    message("The 'counts' assay contains only integers.")
  } else {
    warning("Warning: The 'counts' assay contains decimals. This data is likely already normalized.")
  }
  common_cells <- intersect(colnames(sce), colnames(sce_raw))
  sce_raw <- sce_raw[,common_cells]
  colData(sce_raw) <- colData(sce)[colnames(sce_raw),]
  colData(sce_raw)$celltype = colData(sce_raw)$cell_type
  colData(sce_raw)$batch = colData(sce_raw)$donor_id
  colData(sce_raw)$condition = colData(sce_raw)$disease
  message("celltype:\n", paste(capture.output(table(colData(sce_raw)$celltype)), collapse = "\n"))
  message("batch:\n", paste(capture.output(table(colData(sce_raw)$batch)), collapse = "\n"))
  message("condtion:\n", paste(capture.output(table(colData(sce_raw)$condition)), collapse = "\n"))
  file.remove(temp_h5ad)
  return(sce_raw)
}

				 
human_lung <- function(args){
  getGEOSuppFiles("GSE130148")
  gunzip("GSE130148/GSE130148_raw_counts.RData.gz", remove = FALSE)
  load("GSE130148/GSE130148_raw_counts.RData")
  count_mat = raw_counts
  gunzip("GSE130148/GSE130148_barcodes_cell_types.txt.gz", remove = FALSE)
  coldata =  read.table("GSE130148/GSE130148_barcodes_cell_types.txt", header = TRUE, sep = "\t")
  idx = which(colnames(coldata) == "ID")
  colnames(coldata)[idx] = "batch"
  head(coldata)
  coldata$celltype = as.factor(coldata$celltype)
  coldata$batch = as.factor(coldata$batch)
  message("celltype:\n", paste(capture.output(table(coldata$celltype)), collapse = "\n"))
  message("batch:\n", paste(capture.output(table(coldata$batch)), collapse = "\n"))
  sce = SingleCellExperiment(list(counts = count_mat), colData = coldata)
  return(sce_save)
}

## cellxgene https://cellxgene.cziscience.com/e/9813a1d4-d107-459e-9b2e-7687be935f69.cxg/  https://cellxgene.cziscience.com/collections/5006d6f2-d414-42ed-85d2-d436ee266ac5?explainNewTab
##  download from cellxgene Single-soma transcriptomics of tangle-bearing neurons in Alzheimer’s disease - Inhibitory 
## GSE129308, 6 inhibitory neuron clusters
				 
human_prefrontal_cortex = function(args){
  seurat <- readRDS(url("https://datasets.cellxgene.cziscience.com/4c4dedb7-1f74-4be7-9916-21e00106a1a7.rds"))
  sce = Seurat::as.SingleCellExperiment(seurat)
  colData(sce)$cell = rownames(colData(sce))
  coldata = colData(sce) %>% as_tibble() %>% rename(celltype = Cell.Types) %>% mutate(batch = droplevels(donor_id)) %>%as.data.frame()
  coldata$condition = coldata$disease
  message("celltype:\n", paste(capture.output(table(coldata$celltype)), collapse = "\n"))
  message("batch:\n", paste(capture.output(table(coldata$batch)), collapse = "\n"))
  message("condtion:\n", paste(capture.output(table(coldata$condition)), collapse = "\n"))
  count_mat = counts(sce)[,coldata$cell]
  sce_save = SingleCellExperiment(list(counts = count_mat), colData = coldata)
  return(sce_save)
}			

### human lung atlas https://cellxgene.cziscience.com/e/493a8b60-d676-44d1-b022-d14c1ad0b36c.cxg/
human_lung_atlas = function(args){
  data = readRDS(url("https://datasets.cellxgene.cziscience.com/5311ca08-a915-4bea-a83c-5f2231ba18ef.rds"))
  meta = data@meta.data
  data = GetAssayData(data, assay="RNA", layer="counts")
  data = data[, meta$tissue=="lung"]
  meta = meta[meta$tissue=="lung", ]
  meta$celltype = as.character(droplevels(meta$cell_type))
  meta$batch = as.character(droplevels(meta$donor_id))
  message("celltype:\n", paste(capture.output(table(meta$celltype)), collapse = "\n"))
  message("batch:\n", paste(capture.output(table(meta$batch)), collapse = "\n"))
  meta = meta[,c("celltype", "batch")]
  sce = SingleCellExperiment(assays = list(counts = data),colData = meta)
  return(sce)
}

### macaque retina fovea  GSE118480, from scPSM figshare
macaque_retina_fovea = function(args){
  dir.create("macaque_retina_fovea", recursive = TRUE, showWarnings = FALSE)
  url <- "https://api.figshare.com/v2/file/download/34298483"
  dest <- "macaque_retina_fovea/macaque_retina_fovea.zip"
  
  options(timeout = 3600)
  options(HTTPUserAgent = "Mozilla/5.0")
  
  download.file(
    url = url,
    destfile = dest,
    mode = "wb",
    method = "libcurl"
  )
  #download.file("https://figshare.com/ndownloader/files/34298483","macaque_retina_fovea/macaque_retina_fovea.zip")
  unzip("macaque_retina_fovea/macaque_retina_fovea.zip", exdir="macaque_retina_fovea/")
  retina_expression_matrix <- readRDS("macaque_retina_fovea/retina_expression_matrix.rds")
  retina_metadata <- readRDS("macaque_retina_fovea/retina_metadata.rds")
  
  ## fovea area
  fovea = which(retina_metadata$region == "Fovea")
  coldata_fovea = retina_metadata[fovea,]
  coldata_fovea = coldata_fovea[,-1] %>% mutate(celltype = as.factor(cluster), batch = as.factor(macaque_id) )
  count_fovea = retina_expression_matrix[,rownames(coldata_fovea)]
  message("celltype:\n", paste(capture.output(table(coldata_fovea$celltype)), collapse = "\n"))
  message("batch:\n", paste(capture.output(table(coldata_fovea$batch)), collapse = "\n"))
  sce = SingleCellExperiment(list(counts = count_fovea), colData = coldata_fovea)
  return(sce)
}

## cellxgene https://cellxgene.cziscience.com/e/471647b3-04fe-4c76-8372-3264feb950e8.cxg/
## https://cellxgene.cziscience.com/collections/1cd82f35-026d-48c1-8633-d27ef7485746?explainNewTab
## download from cellxgene CD34+ Fetal Bone Marrow, Fetal Liver, Cord Blood (CITE-seq)
## only use liver
## GSE166895
## Human liver
human_liver = function(args){
  seurat <- readRDS(url("https://datasets.cellxgene.cziscience.com/8e4e400f-30f5-429f-bc38-014e76effe1c.rds"))
  sce = Seurat::as.SingleCellExperiment(seurat)
  colData(sce)$cell = rownames(colData(sce))
  coldata = colData(sce) %>% as_tibble() %>% filter(tissue =="liver") %>% 
    rename(celltype = cell_type) %>% mutate(batch = droplevels(donor_id)) %>%as.data.frame()
  message("celltype:\n", paste(capture.output(table(coldata$celltype)), collapse = "\n"))
  message("batch:\n", paste(capture.output(table(coldata$batch)), collapse = "\n"))
  count_mat = counts(sce)[,coldata$cell]
  sce_save = SingleCellExperiment(list(counts = count_mat), colData = coldata)
  return(sce_save)
}


### mouse Cortex from scPSM figshare (https://figshare.com/articles/dataset/scPSM/19306661)  
###				 (also in https://singlecell.broadinstitute.org/single_cell/study/SCP425/single-cell-comparison-cortex-data)  SCP425
mouse_cortex = function(args){
  dir.create("mouse_cortex", recursive = TRUE, showWarnings = FALSE)
  url <- "https://api.figshare.com/v2/file/download/34292060"
  dest <- "mouse_cortex/mouse_cortex.zip"
  
  options(timeout = 3600)
  options(HTTPUserAgent = "Mozilla/5.0")
  
  download.file(
    url = url,
    destfile = dest,
    mode = "wb",
    method = "libcurl"
  )
  unzip("mouse_cortex/mouse_cortex.zip",exdir="mouse_cortex")
  cortex_expression_matrix <- readRDS("mouse_cortex/cortex_expression_matrix.rds")
  cortex_metadata <- readRDS("mouse_cortex/cortex_metadata.rds")
  ## only use experiment2
  #exp2 = which(cortex_metadata$Experiment=="Cortex2")
  coldata = cortex_metadata# [exp2,]
  no_smart_seq = which(coldata$Method != "Smart-seq2")
  coldata = coldata[no_smart_seq,]
  idx = which(colnames(coldata) == "CellType")
  colnames(coldata)[idx] = "celltype"
  idx = which(colnames(coldata) == "Method")
  colnames(coldata)[idx] = "batch"
  coldata$celltype = as.factor(coldata$celltype)
  coldata$batch = as.factor(coldata$batch)
  message("celltype:\n", paste(capture.output(table(coldata$celltype)), collapse = "\n"))
  message("batch:\n", paste(capture.output(table(coldata$batch)), collapse = "\n"))
  count = cortex_expression_matrix[,rownames(coldata)]
  sce = SingleCellExperiment(list(counts = count),colData = coldata)
  sce 
}	 
				 


				 
				 
				 
