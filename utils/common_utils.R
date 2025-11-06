
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
  if(verbose) message(paste0("Writing: ", file, "."))
  write_h5ad(x.ad, fn, mode = "w", compression = "gzip")
  }else{
    if(verbose) message(paste0("Writing: ", file, "."))
    x$write(
    fn,
    compression='gzip'
    )
  }
  if(verbose) message("Done.")
}

read_normmethod <- function(f) fromJSON(paste(readLines(f), 
                                              collapse=""))$normalize
Full command: /usr/local/lib/R/bin/exec/R --no-save --no-restore --no-echo --no-restore --file=.snakemake/repos/bb2fcfc0ecccc5d1272eb65dfe803f9edc0cdf8125aa21d9a0624395953d892c/run_vizbench.R --args --output_dir out/rawdata/mouse_pancreas/.fc73cbc33384e8af805caeef63502009dc5d075fed448a55aa108e5cf0fc651a/simulate/scdesign3/.9a35770cb788162ac3ef2339501ee66a3442b29594013990097e589e2766f566/normalize/normalize-r/.81dc39bb0dca9db82dd80d2464a44574240337e5dd9d3120a0893477072759db --name mouse_pancreas --simulate.ad out/rawdata/mouse_pancreas/.fc73cbc33384e8af805caeef63502009dc5d075fed448a55aa108e5cf0fc651a/simulate/scdesign3/.9a35770cb788162ac3ef2339501ee66a3442b29594013990097e589e2766f566/mouse_pancreas_simulate.ad --flavour PFlog1pPF --what normalize
Selected category: normalize
Routine selected: PFlog1pPF
Additional parameters: 
name: mouse_pancreas
Verbose: TRUE
location: .snakemake/repos/bb2fcfc0ecccc5d1272eb65dfe803f9edc0cdf8125aa21d9a0624395953d892c
libPaths: /usr/local/lib/R/site-library;/usr/local/lib/R/library
info: sysname=Linux;release=5.4.0-216-generic;version=#236-Ubuntu SMP Fri Apr 11 19:53:21 UTC 2025;nodename=lambda-server;machine=x86_64;login=zhiqian;user=zhiqian;effective_user=zhiqian
Sourcing .. .snakemake/repos/bb2fcfc0ecccc5d1272eb65dfe803f9edc0cdf8125aa21d9a0624395953d892c/utils/common_utils.R
Sourcing .. .snakemake/repos/bb2fcfc0ecccc5d1272eb65dfe803f9edc0cdf8125aa21d9a0624395953d892c/utils/normalize_utils.py
/usr/local/lib/python3.12/dist-packages/scanpy/_utils/__init__.py:35: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('anndata')` instead.
  from anndata import __version__ as anndata_version
/usr/local/lib/python3.12/dist-packages/anndata/__init__.py:70: FutureWarning: Importing read_csv from `anndata` is deprecated. Import anndata.io.read_csv instead.
  return module_get_attr_redirect(attr_name, deprecated_mapping=_DEPRECATED)
/usr/local/lib/python3.12/dist-packages/anndata/__init__.py:70: FutureWarning: Importing read_excel from `anndata` is deprecated. Import anndata.io.read_excel instead.
  return module_get_attr_redirect(attr_name, deprecated_mapping=_DEPRECATED)
/usr/local/lib/python3.12/dist-packages/anndata/__init__.py:70: FutureWarning: Importing read_hdf from `anndata` is deprecated. Import anndata.io.read_hdf instead.
  return module_get_attr_redirect(attr_name, deprecated_mapping=_DEPRECATED)
/usr/local/lib/python3.12/dist-packages/anndata/__init__.py:70: FutureWarning: Importing read_loom from `anndata` is deprecated. Import anndata.io.read_loom instead.
  return module_get_attr_redirect(attr_name, deprecated_mapping=_DEPRECATED)
/usr/local/lib/python3.12/dist-packages/anndata/__init__.py:70: FutureWarning: Importing read_mtx from `anndata` is deprecated. Import anndata.io.read_mtx instead.
  return module_get_attr_redirect(attr_name, deprecated_mapping=_DEPRECATED)
/usr/local/lib/python3.12/dist-packages/anndata/__init__.py:70: FutureWarning: Importing read_text from `anndata` is deprecated. Import anndata.io.read_text instead.
  return module_get_attr_redirect(attr_name, deprecated_mapping=_DEPRECATED)
/usr/local/lib/python3.12/dist-packages/anndata/__init__.py:70: FutureWarning: Importing read_umi_tools from `anndata` is deprecated. Import anndata.io.read_umi_tools instead.
  return module_get_attr_redirect(attr_name, deprecated_mapping=_DEPRECATED)
Sourcing .. .snakemake/repos/bb2fcfc0ecccc5d1272eb65dfe803f9edc0cdf8125aa21d9a0624395953d892c/utils/normalize_utils.R
done running
Converting ReticulateAnnData -> AnnData.Converting AbstractAnnData -> AnnData.Converting R6 -> AnnData.
Writing: out/rawdata/mouse_pancreas/.fc73cbc33384e8af805caeef63502009dc5d075fed448a55aa108e5cf0fc651a/simulate/scdesign3/.9a35770cb788162ac3ef2339501ee66a3442b29594013990097e589e2766f566/normalize/normalize-r/.81dc39bb0dca9db82dd80d2464a44574240337e5dd9d3120a0893477072759db/mouse_pancreas_normalize.ad.
Error in (function (...)  : attempt to apply non-function
Calls: write_ad ... <Anonymous> -> split_named_unnamed -> eval -> eval -> <Anonymous>
Execution halted
