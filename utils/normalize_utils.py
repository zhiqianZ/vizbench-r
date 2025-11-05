import scanpy as sc
import anndata as ad
import numpy as np
import scipy as sp

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
