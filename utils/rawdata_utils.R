
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
  file.remove(temp_h5ad)
  return(sce_raw)
}

human_atheroma <- function(args) {
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
  file.remove(temp_h5ad)
  return(sce_raw)
}


human_glaucoma_pbmc <- function(args) {
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
  file.remove(temp_h5ad)
  return(sce_raw)
}
				 
human_covid_blood <- function(args) {
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
  file.remove(temp_h5ad)
  return(sce_raw)
}


human_colorectal_liver <- function(args) {
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
  file.remove(temp_h5ad)
  return(sce_raw)
}


				 
				 
				 
