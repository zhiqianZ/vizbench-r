
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
  ### single-cell atlas of human pediatric liver (health and IFALD [intestinal failure–associated liver disease])
  url = "https://datasets.cellxgene.cziscience.com/feca90bb-00df-4623-8398-1e3e6a90971d.h5ad"
  temp_h5ad <- tempfile(fileext = ".h5ad")
  options(timeout = 3600)
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
  sce_raw = sce_raw[,sce_raw$celltype != "unkown"]
  file.remove(temp_h5ad)
  return(sce_raw)
}
				 
### huamn pancreas GSM2230757,GSM2230758,GSM2230759,GSM2230760
# file = c("GSM2230757","GSM2230758","GSM2230759","GSM2230760")
# data = lapply(1:4, function(x){
#   getGEOSuppFiles(file[x])
#   curr_data = read.csv(paste0(file[x],"/",file[x],"_human",x,"_umifm_counts.csv.gz"))
#   return(curr_data)
# })
# #data = lapply(1:4, function(i) read.csv(paste0("/home/zhiqian/Benchmark/datasets/HumanPancreas/",file[i])))
# label = lapply(data, function(obj) obj[,3])
# label = unlist(label)
# batch = rep(paste0("inDrop",1:4),unlist(lapply(data,nrow)))
# data = lapply(data,function(obj) {rownames(obj)=obj[,1]; obj[,-(1:3)]})
# count_mat = t(do.call("rbind", data))
# coldata = data.frame(celltype = label, batch = batch)
# rownames(coldata) = colnames(count_mat)
# sce_save = SingleCellExperiment(list(counts = count_mat), colData = coldata)
# unique(coldata$celltype)
# unique(coldata$batch)
# sce_save
# saveRDS(sce_save,"~/Benchmark/Data/sce/human_pancreas.rds")


### Human Lung  GSE130148
# getGEOSuppFiles("GSE130148")
# gunzip("GSE130148/GSE130148_raw_counts.RData.gz", remove = FALSE)
# load("GSE130148/GSE130148_raw_counts.RData")
# count_mat = raw_counts
# gunzip("GSE130148/GSE130148_barcodes_cell_types.txt.gz", remove = FALSE)
# coldata =  read.table("~/Benchmark/raw_data/GSE130148/GSE130148_barcodes_cell_types.txt", header = TRUE, sep = "\t")
# idx = which(colnames(coldata) == "ID")
# colnames(coldata)[idx] = "batch"
# head(coldata)
# coldata$celltype = as.factor(coldata$celltype)
# coldata$batch = as.factor(coldata$batch)
# unique(coldata$celltype)
# unique(coldata$batch)
# sce = SingleCellExperiment(list(counts = count_mat), colData = coldata)
# sce
# saveRDS(sce,"~/Benchmark/Data/sce/human_lung.rds")

### Human PBMC from scPSM figshare (also in https://singlecell.broadinstitute.org/single_cell/study/SCP424/single-cell-comparison-pbmc-data#study-summary)  SCP424
# download.file("https://figshare.com/ndownloader/files/34298480","human_pbmc/human_pbmc.zip")
# unzip("human_pbmc/human_pbmc.zip",exdir="human_pbmc")
# pbmc_expression_matrix <- readRDS("human_pbmc/pbmc_expression_matrix.rds")
# pbmc_metadata = readRDS("human_pbmc/pbmc_metadata.rds")
# ## only use data from experiment 2 to avoid additional layer of batch effect
# exp2 = which(pbmc_metadata$Experiment=="pbmc2")
# coldata = pbmc_metadata[exp2,]
# idx = which(colnames(coldata) == "CellType")
# colnames(coldata)[idx] = "celltype"
# idx = which(colnames(coldata) == "Method")
# colnames(coldata)[idx] = "batch"
# coldata$celltype = as.factor(coldata$celltype)
# coldata$batch = as.factor(coldata$batch)
# unique(coldata$celltype)
# unique(coldata$batch)
# count = pbmc_expression_matrix[,rownames(coldata)]
# sce = SingleCellExperiment(list(counts = count),colData = coldata)
# sce
# saveRDS(sce,"~/Benchmark/Data/sce/human_pbmc.rds")


## cellxgene https://cellxgene.cziscience.com/e/9813a1d4-d107-459e-9b2e-7687be935f69.cxg/  
## https://cellxgene.cziscience.com/collections/5006d6f2-d414-42ed-85d2-d436ee266ac5?explainNewTab
# download from cellxgene Single-soma transcriptomics of tangle-bearing neurons in Alzheimer’s disease - Inhibitory 
# GSE129308, 6 inhibitory neuron clusters
# Human prefrontal cortex
# seurat <- readRDS(url("https://datasets.cellxgene.cziscience.com/4c4dedb7-1f74-4be7-9916-21e00106a1a7.rds"))
# sce = Seurat::as.SingleCellExperiment(seurat)
# colData(sce)$cell = rownames(colData(sce))
# coldata = colData(sce) %>% as_tibble() %>% filter(disease !="Alzheimer disease") %>% rename(celltype = Cell.Types) %>% mutate(batch = droplevels(donor_id)) %>%as.data.frame()
# unique(coldata$celltype)
# unique(coldata$batch)
# count_mat = counts(sce)[,coldata$cell]
# sce_save = SingleCellExperiment(list(counts = count_mat), colData = coldata)
# sce_save
# saveRDS(sce_save,"~/Benchmark/Data/sce/human_prefrontal_cortex.rds")


## cellxgene https://cellxgene.cziscience.com/e/471647b3-04fe-4c76-8372-3264feb950e8.cxg/
# https://cellxgene.cziscience.com/collections/1cd82f35-026d-48c1-8633-d27ef7485746?explainNewTab
# download from cellxgene CD34+ Fetal Bone Marrow, Fetal Liver, Cord Blood (CITE-seq)
# only use liver
# GSE166895
# Human liver
# seurat <- readRDS(url("https://datasets.cellxgene.cziscience.com/8e4e400f-30f5-429f-bc38-014e76effe1c.rds"))
# sce = Seurat::as.SingleCellExperiment(seurat)
# colData(sce)$cell = rownames(colData(sce))
# coldata = colData(sce) %>% as_tibble() %>% filter(tissue =="liver") %>% rename(celltype = cell_type) %>% mutate(batch = droplevels(donor_id)) %>%as.data.frame()
# unique(coldata$celltype)
# unique(coldata$batch)
# count_mat = counts(sce)[,coldata$cell]
# sce_save = SingleCellExperiment(list(counts = count_mat), colData = coldata)
# sce_save
# saveRDS(sce_save,"~/Benchmark/Data/sce/human_liver.rds")

### human lung atlas https://cellxgene.cziscience.com/e/493a8b60-d676-44d1-b022-d14c1ad0b36c.cxg/
# data = readRDS(url("https://datasets.cellxgene.cziscience.com/5311ca08-a915-4bea-a83c-5f2231ba18ef.rds"))
# meta = data@meta.data
# data = GetAssayData(data, assay="RNA", layer="counts")
# data = data[, meta$tissue=="lung"]
# meta = meta[meta$tissue=="lung", ]
# meta$celltype = as.character(droplevels(meta$cell_type))
# meta$batch = as.character(droplevels(meta$donor_id))
# unique(meta$celltype)
# unique(meta$batch)
# meta = meta[,c("celltype", "batch")]
# sce = SingleCellExperiment(assays = list(counts = data),colData = meta)
# sce
# saveRDS(sce,"~/Benchmark/Data/sce/human_lung_atlas.rds")

### macaque retina fovea  GSE118480, from scPSM figshare
# download.file("https://figshare.com/ndownloader/files/34298483","macaque_retina_fovea/macaque_retina_fovea.zip")
# unzip("macaque_retina_fovea/macaque_retina_fovea.zip", exdir="macaque_retina_fovea/")
# retina_expression_matrix <- readRDS("macaque_retina_fovea/retina_expression_matrix.rds")
# retina_metadata <- readRDS("macaque_retina_fovea/retina_metadata.rds")
# 
# ## fovea area
# fovea = which(retina_metadata$region == "Fovea")
# coldata_fovea = retina_metadata[fovea,]
# coldata_fovea = coldata_fovea[,-1] %>% mutate(celltype = as.factor(cluster), batch = as.factor(macaque_id) )
# count_fovea = retina_expression_matrix[,rownames(coldata_fovea)]
# unique(coldata_fovea$celltype)
# length(unique(coldata_fovea$celltype))
# unique(coldata_fovea$batch)
# sce = SingleCellExperiment(list(counts = count_fovea), colData = coldata_fovea)
# sce
# saveRDS(sce,"~/Benchmark/Data/sce/macaque_retina_fovea.rds")
