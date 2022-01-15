import numpy as np
import scipy
import seaborn as sns
import pandas as pd
import genetools


def plot_pca_one_hue_one_timepoint(
    adata,
    ax,
    hue_key,
    title,
    palette,
    continuous_hue=False,
    enable_legends=True,
    marker_size=25,
    all_hue_values_across_timepoints=None,
    annotate_obsnames_on_points=False,
):
    plot_df = adata.obs.assign(
        pca_1=adata.obsm["X_pca"][:, 0], pca_2=adata.obsm["X_pca"][:, 1]
    )
    genetools.plots.scatterplot(
        data=plot_df,
        x_axis_key="pca_1",
        y_axis_key="pca_2",
        hue_key=hue_key,
        continuous_hue=continuous_hue,
        discrete_palette=palette,
        ax=ax,
        marker_size=marker_size,
        alpha=0.9,
        marker="o",
        marker_edge_color="face",  # Enable edges, unless overriden by a HueValueStyle in the palette
        enable_legend=enable_legends,
        legend_hues=all_hue_values_across_timepoints if not continuous_hue else None,
        legend_title=hue_key,
        autoscale=True,
        remove_x_ticks=True,
        remove_y_ticks=True,
        tight_layout=False,
        despine=False,
    )
    ax.set_xlabel("PC 1")
    ax.set_ylabel("PC 2")
    ax.set_title(title)

    # Add annotations - useful for debugging
    if annotate_obsnames_on_points:
        for obsname, position_data in plot_df.iterrows():
            ax.annotate(obsname, (position_data["pca_1"], position_data["pca_2"]))

    return ax


def plot_distance_from_centroid_for_each_hue_value(
    adata_transformed,
    hue_col,
    title,
    ax,
    palette,
    use_rep="X",
):
    # Remove any entries with null values for this hue column
    adata_transformed = adata_transformed[~adata_transformed.obs[hue_col].isna()].copy()

    result = []
    # For each unique hue value in hue_col:
    for hue_value, obs_grp in adata_transformed.obs.groupby(hue_col, observed=True):
        adata_grp = adata_transformed[obs_grp.index]
        # Use raw data or use PCA?
        if use_rep == "X":
            data_to_use = adata_grp.X
        elif use_rep == "X_pca":
            # Use first two dimensions of PCA only, because this is intended to connect to our 2D PCA scatterplots
            data_to_use = adata_grp.obsm["X_pca"][:, :2]
        else:
            raise ValueError("Unrecognized use_rep")

        # Find centroid
        centroid = data_to_use.mean(axis=0)

        # make into 2d array
        centroid = centroid[np.newaxis, :]

        # compute distances from each point to centroid, and flatten into 1d array
        distances = np.ravel(
            scipy.spatial.distance.cdist(data_to_use, centroid, metric="euclidean")
        )

        result.append(pd.Series(distances, index=obs_grp.index))

    # assemble result series with original adata index, but in different order
    series_name = f"Distance from centroid (mean) in {'PC1&2' if use_rep == 'X_pca' else 'raw'} space"
    result = pd.concat(result).rename(series_name)
    adata_transformed.obs = genetools.helpers.merge_into_left(
        adata_transformed.obs, result
    )
    if adata_transformed.obs[series_name].isna().any():
        # sanity checks
        raise ValueError("Distance was not computed for some points")

    sns.boxplot(
        data=adata_transformed.obs,
        x=series_name,
        y=hue_col,
        ax=ax,
        palette=genetools.palette.HueValueStyle.huestyles_to_colors_dict(palette),
    )

    # Add sample size to labels
    ax.set_yticklabels(
        genetools.plots.add_sample_size_to_labels(
            labels=ax.get_yticklabels(), data=adata_transformed.obs, hue_key=hue_col
        )
    )

    # Wrap y-axis text labels
    genetools.plots.wrap_tick_labels(ax, wrap_x_axis=False, wrap_y_axis=True)

    ax.set_title(title)


def make_axis_limits_consistent(axarr):
    """
    make x and y axis limits match across array of axes. set limits to outermost extents encountered in any axes.
    """
    left = min([ax.get_xlim()[0] for ax in axarr])
    right = max([ax.get_xlim()[1] for ax in axarr])

    bottom = min([ax.get_ylim()[0] for ax in axarr])
    top = max([ax.get_ylim()[1] for ax in axarr])

    for ax in axarr:
        ax.set_xlim(left, right)
        ax.set_ylim(bottom, top)

    return axarr, (left, right), (bottom, top)
