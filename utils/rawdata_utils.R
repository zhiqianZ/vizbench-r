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

## -----------------------------------------------------------------------------
## Helper functions
## -----------------------------------------------------------------------------

.warn_if_noninteger_counts <- function(sce) {
  raw_values <- SingleCellExperiment::counts(sce)

  is_integer <- all(raw_values == round(raw_values))

  if (!is_integer) {
    warning(
      "The 'counts' assay contains decimals. ",
      "This data may already be normalized."
    )
  }

  invisible(is_integer)
}

.standardize_coldata <- function(coldata) {
  coldata <- as.data.frame(coldata)

  if (!"celltype" %in% colnames(coldata)) {
    stop("coldata must contain a column named 'celltype'.")
  }

  if (!"batch" %in% colnames(coldata)) {
    stop("coldata must contain a column named 'batch'.")
  }

  coldata$celltype <- as.character(coldata$celltype)
  coldata$batch <- as.character(coldata$batch)

  if ("condition" %in% colnames(coldata)) {
    coldata$condition <- as.character(coldata$condition)
  }

  coldata
}

.finalize_sce <- function(count_mat, coldata) {
  coldata <- .standardize_coldata(coldata)

  if (is.null(rownames(coldata))) {
    if (nrow(coldata) == ncol(count_mat)) {
      rownames(coldata) <- colnames(count_mat)
    } else {
      stop("coldata has no rownames and its nrow does not match ncol(count_mat).")
    }
  }

  common_cells <- intersect(colnames(count_mat), rownames(coldata))

  if (length(common_cells) == 0) {
    if (nrow(coldata) == ncol(count_mat)) {
      rownames(coldata) <- colnames(count_mat)
      common_cells <- colnames(count_mat)
    } else {
      stop("No overlapping cell names between count_mat columns and coldata rownames.")
    }
  }

  count_mat <- count_mat[, common_cells, drop = FALSE]
  coldata <- coldata[common_cells, , drop = FALSE]

  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = count_mat),
    colData = coldata
  )
}

.read_cellxgene_raw <- function(url) {
  options(timeout = 7200)

  temp_h5ad <- tempfile(fileext = ".h5ad")
  on.exit(unlink(temp_h5ad), add = TRUE)

  download.file(url, destfile = temp_h5ad, mode = "wb")

  sce <- zellkonverter::readH5AD(
    file = temp_h5ad,
    X_name = "logcounts",
    raw = TRUE
  )

  sce_raw <- as(
    SingleCellExperiment::altExp(sce, "raw"),
    "SingleCellExperiment"
  )

  assayNames(sce_raw)[assayNames(sce_raw) == "X"] <- "counts"

  .warn_if_noninteger_counts(sce_raw)

  common_cells <- intersect(colnames(sce), colnames(sce_raw))
  sce_raw <- sce_raw[, common_cells]

  colData(sce_raw) <- colData(sce)[colnames(sce_raw), ]

  sce_raw
}

.download_figshare_zip <- function(url, dest, exdir) {
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)

  options(timeout = 3600)
  options(HTTPUserAgent = "Mozilla/5.0")

  download.file(
    url = url,
    destfile = dest,
    mode = "wb",
    method = "libcurl"
  )

  unzip(dest, exdir = exdir)
}

## -----------------------------------------------------------------------------
## Raw data loaders
## -----------------------------------------------------------------------------

mouse_pancreas <- function(args) {
  gsm_ids <- c("GSM2230761", "GSM2230762")

  data_list <- lapply(seq_along(gsm_ids), function(i) {
    tmpdir <- tempdir()

    GEOquery::getGEOSuppFiles(gsm_ids[i], baseDir = tmpdir)

    read.csv(
      file.path(
        tmpdir,
        gsm_ids[i],
        paste0(gsm_ids[i], "_mouse", i, "_umifm_counts.csv.gz")
      )
    )
  })

  celltype <- unlist(lapply(data_list, function(obj) obj[, 3]))
  batch <- rep(
    paste0("inDrop", seq_along(data_list)),
    unlist(lapply(data_list, nrow))
  )

  data_list <- lapply(data_list, function(obj) {
    rownames(obj) <- obj[, 1]
    obj[, -(1:3)]
  })

  count_mat <- t(do.call(rbind, data_list))

  coldata <- data.frame(
    celltype = celltype,
    batch = batch,
    stringsAsFactors = FALSE
  )

  rownames(coldata) <- colnames(count_mat)

  .finalize_sce(count_mat, coldata)
}

human_IFALD_liver <- function(args) {
  url <- "https://datasets.cellxgene.cziscience.com/feca90bb-00df-4623-8398-1e3e6a90971d.h5ad"

  sce_raw <- .read_cellxgene_raw(url)

  colData(sce_raw)$celltype <- colData(sce_raw)$cell_type
  colData(sce_raw)$batch <- colData(sce_raw)$donor_id
  colData(sce_raw)$condition <- colData(sce_raw)$disease

  cond <- as.character(colData(sce_raw)$condition)
  cond[cond != "normal"] <- "intestinal failure-associated liver disease"
  colData(sce_raw)$condition <- cond

  sce_raw <- sce_raw[, sce_raw$celltype != "unknown"]

  sce_raw
}

human_atheroma <- function(args) {
  url <- "https://datasets.cellxgene.cziscience.com/3bab1c0b-d3e3-4a01-840f-d49a8284d989.h5ad"

  sce_raw <- .read_cellxgene_raw(url)

  colData(sce_raw)$celltype <- colData(sce_raw)$cell_type
  colData(sce_raw)$batch <- colData(sce_raw)$donor_id

  sce_raw
}

human_glaucoma_pbmc <- function(args) {
  url <- "https://datasets.cellxgene.cziscience.com/06932380-c04d-4d11-b4af-22de62b031b4.h5ad"

  sce_raw <- .read_cellxgene_raw(url)

  colData(sce_raw)$celltype <- colData(sce_raw)$cell_type
  colData(sce_raw)$batch <- colData(sce_raw)$donor_id
  colData(sce_raw)$condition <- colData(sce_raw)$disease

  sce_raw
}

human_covid_blood <- function(args) {
  url <- "https://datasets.cellxgene.cziscience.com/687c09ff-731a-4e3d-ac07-4c29c33a6338.h5ad"

  sce_raw <- .read_cellxgene_raw(url)

  colData(sce_raw)$celltype <- colData(sce_raw)$cell_type
  colData(sce_raw)$batch <- colData(sce_raw)$donor_id
  colData(sce_raw)$condition <- colData(sce_raw)$disease

  sce_raw
}

human_colorectal_liver <- function(args) {
  url <- "https://datasets.cellxgene.cziscience.com/b425976f-9d73-4388-95dd-e7cd0f8caca0.h5ad"

  sce_raw <- .read_cellxgene_raw(url)

  colData(sce_raw)$celltype <- colData(sce_raw)$cell_type
  colData(sce_raw)$batch <- colData(sce_raw)$donor_id
  colData(sce_raw)$condition <- colData(sce_raw)$disease

  sce_raw
}

human_lung <- function(args) {
  GEOquery::getGEOSuppFiles("GSE130148")

  R.utils::gunzip(
    "GSE130148/GSE130148_raw_counts.RData.gz",
    remove = FALSE,
    overwrite = TRUE
  )

  load("GSE130148/GSE130148_raw_counts.RData")

  count_mat <- raw_counts

  R.utils::gunzip(
    "GSE130148/GSE130148_barcodes_cell_types.txt.gz",
    remove = FALSE,
    overwrite = TRUE
  )

  coldata <- read.table(
    "GSE130148/GSE130148_barcodes_cell_types.txt",
    header = TRUE,
    sep = "\t"
  )

  colnames(coldata)[colnames(coldata) == "ID"] <- "batch"

  coldata$celltype <- as.character(coldata$celltype)
  coldata$batch <- as.character(coldata$batch)

  if (is.null(rownames(coldata)) || !any(rownames(coldata) %in% colnames(count_mat))) {
    if (nrow(coldata) == ncol(count_mat)) {
      rownames(coldata) <- colnames(count_mat)
    }
  }

  .finalize_sce(count_mat, coldata)
}

human_prefrontal_cortex <- function(args) {
  seurat <- readRDS(
    url("https://datasets.cellxgene.cziscience.com/4c4dedb7-1f74-4be7-9916-21e00106a1a7.rds")
  )

  sce <- Seurat::as.SingleCellExperiment(seurat)
  colData(sce)$cell <- rownames(colData(sce))

  coldata <- colData(sce) |>
    as.data.frame() |>
    dplyr::rename(celltype = Cell.Types) |>
    dplyr::mutate(
      batch = as.character(droplevels(donor_id)),
      condition = as.character(disease)
    )

  rownames(coldata) <- coldata$cell

  count_mat <- SingleCellExperiment::counts(sce)[, rownames(coldata), drop = FALSE]

  .finalize_sce(count_mat, coldata)
}

human_lung_atlas <- function(args) {
  seurat <- readRDS(
    url("https://datasets.cellxgene.cziscience.com/5311ca08-a915-4bea-a83c-5f2231ba18ef.rds")
  )

  meta <- seurat@meta.data
  count_mat <- Seurat::GetAssayData(seurat, assay = "RNA", layer = "counts")

  keep <- meta$tissue == "lung"

  meta <- meta[keep, , drop = FALSE]
  count_mat <- count_mat[, rownames(meta), drop = FALSE]

  coldata <- data.frame(
    celltype = as.character(droplevels(meta$cell_type)),
    batch = as.character(droplevels(meta$donor_id)),
    stringsAsFactors = FALSE
  )

  rownames(coldata) <- rownames(meta)

  .finalize_sce(count_mat, coldata)
}

macaque_retina_fovea <- function(args) {
  dir_name <- "macaque_retina_fovea"
  zip_path <- file.path(dir_name, "macaque_retina_fovea.zip")

  .download_figshare_zip(
    url = "https://api.figshare.com/v2/file/download/34298483",
    dest = zip_path,
    exdir = dir_name
  )

  retina_expression_matrix <- readRDS(
    file.path(dir_name, "retina_expression_matrix.rds")
  )

  retina_metadata <- readRDS(
    file.path(dir_name, "retina_metadata.rds")
  )

  keep <- retina_metadata$region == "Fovea"

  coldata <- retina_metadata[keep, , drop = FALSE]
  coldata <- coldata[, -1, drop = FALSE] |>
    dplyr::mutate(
      celltype = as.character(cluster),
      batch = as.character(macaque_id)
    )

  count_mat <- retina_expression_matrix[, rownames(coldata), drop = FALSE]

  .finalize_sce(count_mat, coldata)
}

human_liver <- function(args) {
  seurat <- readRDS(
    url("https://datasets.cellxgene.cziscience.com/8e4e400f-30f5-429f-bc38-014e76effe1c.rds")
  )

  sce <- Seurat::as.SingleCellExperiment(seurat)
  colData(sce)$cell <- rownames(colData(sce))

  coldata <- colData(sce) |>
    as.data.frame() |>
    dplyr::filter(tissue == "liver") |>
    dplyr::rename(celltype = cell_type) |>
    dplyr::mutate(
      batch = as.character(droplevels(donor_id))
    )

  rownames(coldata) <- coldata$cell

  count_mat <- SingleCellExperiment::counts(sce)[, rownames(coldata), drop = FALSE]

  .finalize_sce(count_mat, coldata)
}

mouse_cortex <- function(args) {
  dir_name <- "mouse_cortex"
  zip_path <- file.path(dir_name, "mouse_cortex.zip")

  .download_figshare_zip(
    url = "https://api.figshare.com/v2/file/download/34292060",
    dest = zip_path,
    exdir = dir_name
  )

  cortex_expression_matrix <- readRDS(
    file.path(dir_name, "cortex_expression_matrix.rds")
  )

  cortex_metadata <- readRDS(
    file.path(dir_name, "cortex_metadata.rds")
  )

  coldata <- cortex_metadata
  coldata <- coldata[coldata$Method != "Smart-seq2", , drop = FALSE]

  colnames(coldata)[colnames(coldata) == "CellType"] <- "celltype"
  colnames(coldata)[colnames(coldata) == "Method"] <- "batch"

  coldata <- coldata[coldata$celltype != "Unassigned", , drop = FALSE]

  coldata$celltype <- as.character(coldata$celltype)
  coldata$batch <- as.character(coldata$batch)

  count_mat <- cortex_expression_matrix[, rownames(coldata), drop = FALSE]

  .finalize_sce(count_mat, coldata)
}
