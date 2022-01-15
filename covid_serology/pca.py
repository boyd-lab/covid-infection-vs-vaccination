import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import genetools
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, FunctionTransformer
from sklearn.decomposition import PCA
from . import config


def _fit_transformations(adata, random_state=0, svd_solver="arpack", **kwargs):
    """Fit pre-processing and PCA on an adata from one time point.

    This transformation should approximately replicate:

        sc.pp.log1p(adata)
        sc.pp.scale(adata)
        sc.tl.pca(adata, n_comps=n_pcs)

    Except it does not modify adata in place.
    Instead this returns a transformer that can be used to modify any anndata.
    """
    # log1p, scale columns, and run PCA
    pipeline = make_pipeline(
        FunctionTransformer(np.log1p, validate=True),
        StandardScaler(),
        PCA(
            random_state=random_state,
            svd_solver=svd_solver,
            **kwargs,
        ),
    )
    return pipeline.fit(adata.X)


def _apply_transformations(adata, pipeline):
    """Apply pre-fit pre-processing and PCA to an adata from a new time point."""
    # run log1p and scaler only
    X_log1p = pipeline.named_steps["functiontransformer"].transform(adata.X)
    X_scaled = pipeline.named_steps["standardscaler"].transform(X_log1p)

    # run log1p, scaler, and PCA all together
    X_pca = pipeline.transform(adata.X)

    ## update adata object:
    # store original data in raw
    adata.raw = adata

    # don't store log1p data anywhere

    # store scaled in X
    adata.X = X_scaled

    # store PCA info
    adata.obsm["X_pca"] = X_pca
    pca_ = pipeline.named_steps["pca"]
    adata.uns["pca"] = {
        "variance": pca_.explained_variance_,
        "variance_ratio": pca_.explained_variance_ratio_,
        "params": {"zero_center": True, "use_highly_variable": False},
    }
    adata.varm["PCs"] = pca_.components_.T

    return adata


def cleanup_adata(adata):
    """After filtering an anndata, remove missing entries. Returns a copy."""
    adata = adata.copy()

    # remove measurements (columns) that aren't defined for any subject (i.e. blank in all rows)
    # in scanpy lingo, cells are rows and genes are columns, meaning filter_genes removes columns
    sc.pp.filter_genes(adata, min_cells=1)

    # remove subjects (rows) with incomplete entries for this timepoint, i.e. remove any rows that have any NaNs in the remaining columns
    sc.pp.filter_cells(adata, min_genes=adata.shape[1])

    # remove measurements (columns) that aren't defined for all remaining subjects, i.e. remove any columns that have any NaNs among the remaining rows.
    sc.pp.filter_genes(adata, min_cells=adata.shape[0])

    return adata


def export_pca_diagnostics(adata, output_prefix):
    adata.var_names.to_series().to_csv(
        f"{config.paths.output_dir}/{output_prefix}.var_names.txt", index=None
    )

    # how many PCs do we want
    sc.pl.pca_variance_ratio(adata, log=True, show=False)
    fig = plt.gcf()
    genetools.plots.savefig(
        fig, f"{config.paths.output_dir}/{output_prefix}.pca_variance_ratio.png", dpi=72
    )
    plt.close(fig)

    # explained variance ratio
    explained_variance_ratio_df = pd.DataFrame(
        [
            {"number": f"PC {i+1}", "variance explained": f"{round(v * 100, 2)}%"}
            for i, v in enumerate(adata.uns["pca"]["variance_ratio"])
        ]
    ).set_index("number")
    explained_variance_ratio_df.to_csv(
        f"{config.paths.output_dir}/{output_prefix}.pc_variance_explained.tsv", sep="\t"
    )

    # PC1, PC2 loadings
    pc_loadings = (
        pd.DataFrame(adata.varm["PCs"], index=adata.var_names)[[0, 1]]
        .rename(columns={0: "loading_pc1", 1: "loading_pc2"})
        .sort_values("loading_pc1", ascending=False)
    )
    pc_loadings.to_csv(
        f"{config.paths.output_dir}/{output_prefix}.pc_loadings.tsv", sep="\t"
    )


def pca_one_timepoint_against_reference(adata, output_prefix, transformer=None):
    """PCA this timepoint using PCA transformation fit on another timepoint (optional)."""
    if transformer is None:
        transformer = _fit_transformations(adata)
    adata = _apply_transformations(adata, transformer)
    export_pca_diagnostics(adata, output_prefix)
    return adata
