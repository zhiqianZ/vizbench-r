load_pkgs <- function() {
  library(scDesign3)
  library(SingleCellExperiment)
  library(Seurat)
  library(pbmcapply)
  library(parallel)
}

## NOTE: code below uses 40 cores!!

# later (at metric calculation, we'll need some parameters)
# instead of recomputing, this should be calculated at simulation time and read here.
# for(data in Dataset){
#     simulation = readRDS(file.path(data_dir,"Data/Simulation",
#                                    data,paste0(dataset_name, "_simulation.rds")))
#     para = readRDS(file.path(data_dir,"Data/Simulation",
#                              data,paste0(dataset_name, "_para.rds")))
#     celltype = simulation$meta$celltype
#     batch = simulation$meta$batch
#     ls = log10(colSums(simulation$counts))
#     zp <- colMeans(simulation$counts==0)
#     rm(simulation)
#     idx = lapply(unique(batch), function(b) tapply(which(batch==b), celltype[batch==b], function(id) id[1]))
#     names(idx) = unique(batch)
#     mean_par = lapply(idx, function(idx) para$mean_mat[idx,])
#     mean_par = lapply(mean_par, function(obj) {rownames(obj) = names(idx[[1]]); obj})
#     
#     variance_par = lapply(idx, function(idx)  para$mean_mat[idx,]^2 * para$sigma_mat[idx,])
#     variance_par = lapply(variance_par, function(obj) {rownames(obj) = names(idx[[1]]); obj})
#     
#     rm(para)
#     for(norm in Normalization){
#       for(integ in Integration){
#         for(visual in Visualization){
#           res = try(readRDS(paste0(path_visualization_results,"/",data,"_",norm,"+",integ,"+",visual,".rds")))
#           for(met in Metric){
#             if(!inherits(res, "try-error")){
#               val = Metric_eval(data =res, 
#                                 mean_par = mean_par, variance_par = variance_par,
#                                 library_size = ls, zero_proportion = zp,
#                                 celltype = celltype, batch = batch,
#                                 metric = met, n.cores = 20) 
#               metric_df[k, ] = c(data, norm, integ, visual, met, val)
#             }else{
#               metric_df[k, ] = c(data, norm, integ, visual, met, NA)
#             }
#             print(k)
#             saveRDS(metric_df,paste0(path_save_evaluation,"/",data,"_evaluation_metrics.rds"))
#             k=k+1
#           }
#         }
#       }
#     }
#   }


scdesign3 <- function(args) {
  
  sce <- read_sce(args$rawdata.ad)
  data = counts(sce)
  coldat = colData(sce)
  
  # QC
  batch = unique(coldat$batch)
  meta = data.frame(cbind(celltype=as.character(coldat$celltype),
                          batch=as.character(coldat$batch)))
  batch_filtered = colnames(table(meta))[which(table(meta)<5&table(meta)>0, arr.ind = T)[,2]]
  celltype_filtered = rownames(table(meta))[which(table(meta)<5&table(meta)>0, arr.ind = T)[,1]]
  filtered_genes = which(apply(data,1,function(x) length(which(x!=0)))<10)
  if(length(filtered_genes)>0){
    data = data[-filtered_genes,]  
  }
  if(length(batch_filtered) > 0){
    filter = paste0(batch_filtered,"_", celltype_filtered)
    idx = which(paste0(meta$batch,"_", meta$celltype) %in% filter)
    data = data[, -idx]
    meta = meta[-idx, ]
    coldat = coldat[-idx,]
  }
  
  # HVG 
  seurat = CreateSeuratObject(data)
  seurat = FindVariableFeatures(seurat,nfeatures = 2000) 
  # MR: 1000 to make it faster?
  hvg = VariableFeatures(seurat)
  
  data = data[hvg,]
  if(args$verbose) message(paste0(dim(data), collapse=","))
  ## scDesign3
  sce <- SingleCellExperiment(assay = list(counts = data), 
                              colData = coldat)

  set.seed(123)
  data <- construct_data(
    sce = sce,
    assay_use = "counts",
    celltype = "celltype",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = "batch",
    corr_by = "celltype",
    ncell = 3000
    #ncell = 100000
  )
  
  if(args$verbose) message(mean(data$count_mat>0))
  if(args$verbose) message(paste0(dim(data$count_mat), collapse=","))
  
  marginal <- fit_marginal(
    data = data,
    predictor = "gene",
    #mu_formula = "celltype+s(batch,bs='re')",
    mu_formula = "celltype+batch",
    sigma_formula = "celltype",
    family_use = "nb",
    n_cores = 10,
    usebam = FALSE,
    parallelization = "mcmapply",
    # parallelization = "pbmcmapply",
    trace = TRUE
  )
  
  copula <- fit_copula(
    sce = sce,
    assay_use = "counts",
    marginal_list = marginal,
    family_use = "nb",
    copula = "gaussian",
    n_cores = 10,
    input_data = data$dat
  )
  
  para <- extract_para(
    sce = sce,
    marginal_list = marginal,
    n_cores = 10,
    family_use = "nb",
    new_covariate = data$newCovariate,
    parallelization = "mcmapply",
    data = data$dat#,
    #parallelization = "pbmcmapply"
  )
  
  # should we save the parameters? atm, just saving the dataset
  # dir.create(file.path(base_dir, "Data/Simulation"))
  # dir.create(file.path(base_dir, "Data/Simulation/Pilot"))
  # saveRDS(para,file.path(base_dir, "Data/Simulation/Pilot/MousePancreas_para.rds"))

  if(args$verbose) message(paste0(names(data), collapse=","))
  
  set.seed(123)
  newcounts = simu_new(
    sce = sce,
    filtered_gene = data$filtered_gene,
    mean_mat = para$mean_mat,
    sigma_mat = para$sigma_mat,
    zero_mat = para$zero_mat,
    quantile_mat = NULL,
    copula_list = copula$copula_list,
    n_cores = 10,
    family_use = "nb",
    input_data = data$dat,
    new_covariate = data$newCovariate,
    important_feature = copula$important_feature,
    parallelization = "pbmcmapply"
  )

  seurat.obj = CreateSeuratObject(newcounts, 
				  meta.data=data$newCovariate)
  return(seurat.obj)
}

