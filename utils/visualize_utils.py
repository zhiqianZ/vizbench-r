import scanpy as sc
import anndata as ad
import numpy as np
import scipy as sp
import os


def _neighbors_umap(latent, n_neighbors, min_dist, n_jobs=1, seed=42):
    """Core: kNN graph + UMAP layout on a raw cells x dims matrix."""
    latent = np.asarray(latent, dtype=np.float64)

    adata = ad.AnnData(
        X=np.zeros((latent.shape[0], 1), dtype=np.float32),  # dummy X
        obsm={"integrated": latent},
    )

    sc.pp.neighbors(
        adata,
        use_rep="integrated",
        n_neighbors=int(n_neighbors),
        random_state=int(seed),
    )
    sc.tl.umap(
        adata,
        min_dist=float(min_dist),
        random_state=int(seed),
    )
    return np.asarray(adata.obsm["X_umap"])


def scanpyUMAP(args):
    """Baseline scanpy UMAP. Unoptimized -- scanpy defaults."""
    print("Running scanpyUMAP")
    adata = sc.read_h5ad(args["integrate.ad"])
    npcs = int(args["npcs"])

    latent = np.asarray(adata.obsm["integrated"])[:, :npcs]   # explicit slice
    vis = _neighbors_umap(
        latent,
        n_neighbors=15,    # scanpy default
        min_dist=0.5,      # scanpy default
        n_jobs=int(args.get("nthreads", 1)),
    )
    return vis
    
def scanpyUMAP(args):
    """Baseline scanpy UMAP. Unoptimized -- scanpy defaults."""
    print("Running scanpyUMAP")

    adata = sc.read_h5ad(args["integrate.ad"])
    npcs = int(args["npcs"])

    print("adata.n_obs:", adata.n_obs)
    print("integrated shape:", adata.obsm["integrated"].shape)

    latent = np.asarray(
        adata.obsm["integrated"]
    )[:, :npcs]

    print("latent shape:", latent.shape)
    print("non-finite rows:",
          np.sum(~np.isfinite(latent).all(axis=1)))

    assert latent.shape[0] == adata.n_obs

    vis = _neighbors_umap(
        latent,
        n_neighbors=15,
        min_dist=0.5,
        n_jobs=int(args.get("nthreads", 1)),
    )

    vis = np.asarray(vis)

    print("vis shape:", vis.shape)

    assert vis.ndim == 2
    assert vis.shape[0] == adata.n_obs, (
        f"UMAP returned {vis.shape[0]} rows, "
        f"expected {adata.n_obs}"
    )

    return vis


def scanpy_umap_from_matrix(latent, n_neighbors, min_dist, n_jobs=1, seed=42):
    """
    Optimizer adapter. `latent` is already sliced to npcs columns by the caller.
    Called once per (original, permuted) x grid-row from the R engine.
    """
    return _neighbors_umap(
        latent,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_jobs=n_jobs,
        seed=seed,
    )


def scanpy_umap_grid(latent, grid, n_jobs=1, seed=42):
    """
    Optional: run a whole grid in one call, caching the kNN graph.

    sc.pp.neighbors depends only on n_neighbors; sc.tl.umap depends on min_dist.
    So we build the graph once per n_neighbors and lay it out once per min_dist,
    instead of rebuilding the graph for every (n_neighbors, min_dist) pair.

    `grid` is a list of (n_neighbors, min_dist) tuples.
    Returns a dict keyed by (n_neighbors, min_dist) -> cells x 2 array.
    """
    latent = np.asarray(latent, dtype=np.float64)

    by_k = {}
    for k, m in grid:
        by_k.setdefault(int(k), []).append(float(m))

    out = {}
    for k, min_dists in by_k.items():
        adata = ad.AnnData(
            X=np.zeros((latent.shape[0], 1), dtype=np.float32),
            obsm={"integrated": latent},
        )
        sc.pp.neighbors(
            adata, use_rep="integrated",
            n_neighbors=k, random_state=int(seed),
        )
        for m in min_dists:                      # graph reused across min_dist
            sc.tl.umap(adata, min_dist=m, random_state=int(seed))
            out[(k, m)] = np.asarray(adata.obsm["X_umap"]).copy()

    return out

"""
def scanpyUMAP(args):
  print("Running scanpyUMAP")
  adata_path = args["integrate.ad"]
  npcs = args['npcs']
  adata = sc.read_h5ad(adata_path)
  sc.pp.neighbors(adata, n_pcs=npcs, use_rep = "integrated")
  sc.tl.umap(adata)
  return adata.obsm['X_umap']
  
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
