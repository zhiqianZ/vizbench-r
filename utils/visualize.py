import scanpy as sc
import anndata as ad
import numpy as np
import scipy as sp

def scanpyUMAP(args):
  print("Running scanpyUMAP")
  adata_path = args["integrate.ad"]
  adata = sc.read_h5ad(adata_path)
  sc.pp.neighbors(adata, n_pcs=50, use_rep = "integrated")
  sc.tl.umap(adata)
  adata.obsm['X_umap']
  
    
def graphFA(args):
  print("Running graphFA")
  adata_path = args["integrate.ad"]
  adata = sc.read_h5ad(adata_path)
  sc.pp.neighbors(adata, n_pcs=50, use_rep = "integrated")
  sc.tl.draw_graph(adata)
  adata.obsm['X_draw_graph_fa']
  
