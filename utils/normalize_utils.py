import scanpy as sc
import anndata as ad
import numpy as np
import scipy as sp

def log1pCPMedian(args):
    print("Running log1pCPMedian")
    adata_path = args["simulate.ad"]
    adata = sc.read_h5ad(adata_path)
    
    X_data = sp.sparse.lil_matrix(adata.shape, dtype=np.float32)
    # Batch-wise normalization
    for batch in adata.obs["batch"].unique():
        print(f"Normalizing batch: {batch}")
        batch_idx = adata.obs["batch"] == batch
        idx = np.where(batch_idx)[0]
        # Create lightweight AnnData for transformation
        sub = sc.AnnData(X=adata.layers['counts'][idx, :])
        sc.pp.normalize_total(sub, target_sum=1e4, inplace=True)
        sc.pp.log1p(sub)
        for i, orig_idx in enumerate(idx):
            X_data[orig_idx, :] = sub.X[i, :]
    
    # Store in compressed format
    adata.layers["data"] = X_data.tocsr()
    
    return adata


def log1pPF(args):
    print("Running log1pPF")

    adata_path = args["simulate.ad"]
    adata = sc.read_h5ad(adata_path)

    n_cells, n_genes = adata.shape
    X_data = np.zeros((n_cells, n_genes), dtype=np.float32)

    for batch in adata.obs["batch"].unique():
        print(f"Normalizing batch: {batch}")
        batch_mask = adata.obs["batch"] == batch
        idx = np.where(batch_mask)[0]

        # Subset just counts layer and normalize
        X_data[idx, :] = logPF(adata.layers["counts"][idx, :]).toarray()
        
    adata.layers["data"] = sp.sparse.csr_matrix(X_data)

    return adata

def PFlog1pPF(args):
    print("Running PFlog1pPF")

    adata_path = args["simulate.ad"]
    adata = sc.read_h5ad(adata_path)

    n_cells, n_genes = adata.shape
    X_data = np.zeros((n_cells, n_genes), dtype=np.float32)

    for batch in adata.obs["batch"].unique():
        print(f"Normalizing batch: {batch}")
        batch_mask = adata.obs["batch"] == batch
        idx = np.where(batch_mask)[0]

        # Subset just counts layer and normalize
        X_data[idx, :] = PFlogPF(adata.layers["counts"][idx, :]).toarray()
        
    adata.layers["data"] = sp.sparse.csr_matrix(X_data)

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
  return norm
 
def PFlogPF(counts):
  norm = norm_pf_log_pf(sp.sparse.csr_matrix(counts))
  return norm
