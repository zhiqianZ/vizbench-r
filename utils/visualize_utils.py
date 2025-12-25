import scanpy as sc
import anndata as ad
import numpy as np
import scipy as sp
import os

def scanpyUMAP(args):
  print("Running scanpyUMAP")
  adata_path = args["integrate.ad"]
  npcs = args['npcs']
  adata = sc.read_h5ad(adata_path)
  sc.pp.neighbors(adata, n_pcs=npcs, use_rep = "integrated")
  sc.tl.umap(adata)
  return adata.obsm['X_umap']
  
"""
def graphFA(args):
  print("Running graphFA")
  print("== Environment thread variables ==")
  for k in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "OMP_DYNAMIC",
    "MKL_DYNAMIC",
  ]:
    print(f"{k} = {os.getenv(k)}")
  adata_path = args["integrate.ad"]
  npcs = args['npcs']
  adata = sc.read_h5ad(adata_path)
  sc.pp.neighbors(adata, n_pcs=npcs, use_rep = "integrated")
  sc.tl.draw_graph(adata)
  return adata.obsm['X_draw_graph_fa']
  """
