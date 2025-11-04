import scanpy as sc
import anndata as ad
import numpy as np
import scipy as sp

def log1pCPMedian(args):
    print("Running log1pCPMedian")

    adata_path = args['simulate_ad'] if isinstance(args, dict) else args.simulate_ad
    adata = sc.read_h5ad(adata_path)

    X_norm = np.zeros_like(adata.X, dtype=np.float32)
    for batch in adata.obs['batch'].unique():
        print(f"Normalizing batch: {batch}")
        batch_idx = adata.obs['batch'] == batch
        sub = adata[batch_idx].copy()
        sc.pp.normalize_total(sub, target_sum=None, inplace=True)
        sc.pp.log1p(sub)
        X_norm[batch_idx, :] = sub.X
    adata.layers["data"] = X_norm
  
    return adata

def log1pPF(args):
    print("Running log1pPF")
    adata_path = args['simulate_ad'] if isinstance(args, dict) else args.simulate_ad
    adata = ad.read_h5ad(adata_path)

    count_layers = [k for k in adata.layers.keys() if 'counts' in k.lower()]
    for layer in count_layers:
        new_layer = layer.replace('counts', 'data')
        X = adata.layers[layer]
        norm = logPF(X)
        adata.layers[new_layer] = norm
      
    return adata


def PFlog1pPF(adata):
    print("Running PFlog1pPF")

    count_layers = [k for k in adata.layers.keys() if 'counts' in k.lower()]
    for layer in count_layers:
        new_layer = layer.replace('counts', 'data')
        X = adata.layers[layer]
        norm = PFlogPF(X)
        adata.layers[new_layer] = norm
      
    return adata

#### Contribute to xxx
def do_pf(mtx, sf = None):
    pf = mtx.sum(axis=1).A.ravel()
    if not sf:
        sf = pf.mean()
    pf = sp.sparse.diags(sf/pf) @ mtx
    return pf

def norm_pf(mtx):
    return do_pf(mtx)

def norm_log_pf(mtx):
    pf_log = np.log1p(do_pf(mtx))
    return pf_log

def norm_pf_log_pf(mtx):
    pf_log_pf = do_pf(np.log1p(do_pf(mtx)))
    return pf_log_pf
   
def logPF(counts):
  norm = norm_log_pf(sp.sparse.csr_matrix(counts))
  return norm.T
 
def PFlogPF(counts):
  norm = norm_pf_log_pf(sp.sparse.csr_matrix(counts))
  return norm.T
