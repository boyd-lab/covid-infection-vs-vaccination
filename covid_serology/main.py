import numpy as np
import matplotlib.pyplot as plt
import genetools
import sys
from slugify import slugify
from . import config, pca, plots


def transform_and_plot_pcas_multiple_timepoints_against_common_reference(
    adata,
    plot_name,
    reference_timepoint,
    colors_to_plot,
    palette_dict,
    enforce_cutoffs=False,
    only_timepoints_with_multiple_source_cohorts=True,
    timepoint_sort_order=None,
    annotate_obsnames_on_points=False,
):
    """
    For each timepoint: transform against reference timepoint's PCA (computed here), and plot.

    if enforce_cutoffs is True: at each timepoint, only choose subjects who don't violate any plate-11 cutoffs at that timepoint.
    otherwise include all rows.
    """
    timepoint_adatas_pcaed = transform_pcas_multiple_timepoints_against_common_reference(
        adata=adata,
        plot_name=plot_name,
        reference_timepoint=reference_timepoint,
        enforce_cutoffs=enforce_cutoffs,
        timepoint_sort_order=timepoint_sort_order,
        only_timepoints_with_multiple_source_cohorts=only_timepoints_with_multiple_source_cohorts,
    )
    plot_pcas_multiple_timepoints_against_common_reference(
        timepoint_adatas_pcaed=timepoint_adatas_pcaed,
        colors_to_plot=colors_to_plot,
        reference_timepoint=reference_timepoint,
        plot_name=plot_name,
        palette_dict=palette_dict,
        annotate_obsnames_on_points=annotate_obsnames_on_points,
    )
    return timepoint_adatas_pcaed


def transform_pcas_multiple_timepoints_against_common_reference(
    adata,
    plot_name,
    reference_timepoint,
    enforce_cutoffs=False,
    timepoint_sort_order=None,
    only_timepoints_with_multiple_source_cohorts=True,
):
    """
    Run PCA on a time course, with one timepoint used to train the PCA, and all other timepoints projected onto that PCA reference:

    First, fit a transformation on the reference timepoint:
    - Filter the Subjects x Measurements matrix to measurements from that timepoint
    - If we are enforcing cutoffs: only choose subjects who don't violate any p11 cutoffs at that timepoint. Otherwise use all rows.
    - Remove missing entries from the matrix after filtering (e.g. columns that are now blank)
    - Create a pipeline of `log1p`, `StandardScaler`, and `PCA`.
    - Fit that pipeline on the samples x measurements matrix.

    Then apply that PCA transformation to each timepoint, i.e. at each timepoint:
    - Filter the Subjects x Measurements matrix to measurements from that timepoint
    - If we are enforcing cutoffs: only choose subjects who don't violate any p11 cutoffs at that timepoint. Otherwise use all rows.
    - Remove missing entries from the matrix after filtering (e.g. columns that are now blank)
    - Skip any timepoints with fewer than 2 source cohorts remaining or with an empty matrix after the cleanup
    - Run the stored `log1p`, `StandardScaler`, `PCA` transformations

    `enforce_cutoffs` arugment: at each timepoint, only choose subjects who don't violate any p11 cutoffs at that timepoint. otherwise choose all rows.
    """

    # create PCA reference on reference timepoint
    # note: slice(None) is the same as : slicer
    adata_reference_subset = adata[
        (adata.obs[f"any_p11_cutoffs_violated:{reference_timepoint}"] == False)
        if enforce_cutoffs
        else slice(None),
        adata.var["timepoint"] == reference_timepoint,
    ]
    transformer = pca._fit_transformations(
        pca.cleanup_adata(adata_reference_subset),
    )

    # transform all timepoints using that PCA reference
    timepoint_adatas_pcaed = []
    for timepoint_ix, (timepoint, var_by_timepoint) in enumerate(
        adata.var.groupby("timepoint")
    ):
        row_slice = slice(None)
        if enforce_cutoffs:
            row_slice = adata.obs[f"any_p11_cutoffs_violated:{timepoint}"] == False
            sys.stderr.write(
                f"Timepoint {timepoint}: removing {adata.shape[0] - adata[row_slice].shape[0]} rows due to plate 11 cutoffs.\n"
            )
        adata_subset = adata[row_slice, var_by_timepoint.index]

        adata_subset = pca.cleanup_adata(adata_subset)
        if adata_subset.shape[0] == 0 or adata_subset.shape[1] == 0:
            # adata is now empty. don't proceed
            sys.stderr.write(f"Empty dataset for timepoint {timepoint} - skipping.\n")
            continue
        if (
            only_timepoints_with_multiple_source_cohorts
            and adata_subset.obs["source_cohort"].nunique() < 2
        ):
            sys.stderr.write(
                f"Fewer than two source cohorts for timepoint {timepoint} - skipping.\n"
            )
            continue

        # convert plot name to valid filename
        output_prefix = slugify(
            f"{plot_name}.{timepoint}.on_common_reference.{reference_timepoint}"
        )

        # Apply log-transform + scaling + PCA that was fit on reference data previously
        adata_transformed = pca.pca_one_timepoint_against_reference(
            adata_subset,
            output_prefix=output_prefix,
            transformer=transformer,
        )

        timepoint_adatas_pcaed.append((timepoint, adata_transformed))

    if timepoint_sort_order is not None:
        observed_timepoints = [timepoint for (timepoint, _) in timepoint_adatas_pcaed]
        # filter sort order down to timepoints we actually have
        timepoint_sort_order = [
            timepoint
            for timepoint in timepoint_sort_order
            if timepoint in observed_timepoints
        ]
        # sort timepoint_adatas_pcaed according to desired sort order
        # note: can't unpack tuples in lambda declaration: https://stackoverflow.com/q/21892989/130164 and https://www.python.org/dev/peps/pep-3113/
        timepoint_adatas_pcaed.sort(
            key=lambda timepoint_and_adata_tuple: timepoint_sort_order.index(
                timepoint_and_adata_tuple[0]
            )
        )

    return timepoint_adatas_pcaed


def plot_pcas_multiple_timepoints_against_common_reference(
    timepoint_adatas_pcaed,
    colors_to_plot,
    reference_timepoint,
    plot_name,
    palette_dict,
    annotate_obsnames_on_points=False,
):
    """
    Plot one figure per hue of interest:
    PC1 vs PC2 for each timepoint on a consistent PCA reference,
    with each timepoint on a separate subplot (with consistent axes limits).

    Also plot distribution of distances from each hue group to its centroid.
    """
    n_timepoints = len(timepoint_adatas_pcaed)

    # Plot each hue on PCA for each timepoint
    for (hue_col, is_continuous) in colors_to_plot:

        all_hue_values_across_timepoints = None
        if not is_continuous:
            # to construct legends, get all values this discrete hue column takes across timepoints
            # scan across all timepoint anndatas
            all_hue_values_across_timepoints = set(
                np.hstack(
                    [
                        ad.obs[hue_col].unique()
                        for (timepoint, ad) in timepoint_adatas_pcaed
                    ]
                )
            )

        # set up time course figure for this hue of interest. the figure is a single row of subplots (dim: # of timepoints)
        fig, axarr = plt.subplots(nrows=1, ncols=n_timepoints, figsize=(32, 5))

        # plot for each timepoint
        for timepoint_ix, ((timepoint, adata_transformed), ax) in enumerate(
            zip(timepoint_adatas_pcaed, axarr)
        ):
            plots.plot_pca_one_hue_one_timepoint(
                adata=adata_transformed,
                ax=ax,
                hue_key=hue_col,
                title=(
                    f"{timepoint} (reference)"
                    if timepoint == reference_timepoint
                    else timepoint
                ),
                enable_legends=(timepoint_ix == n_timepoints - 1),
                palette=palette_dict,
                continuous_hue=is_continuous,
                marker_size=25,
                all_hue_values_across_timepoints=all_hue_values_across_timepoints,
                annotate_obsnames_on_points=annotate_obsnames_on_points,
            )

        # make all axis limits consistent
        plots.make_axis_limits_consistent(axarr)

        # save time course figure
        # convert plot name to valid filename
        plot_fname = (
            f"{config.paths.output_dir}/{slugify(plot_name)}.pca.{slugify(hue_col)}"
        )
        # rasterized
        genetools.plots.savefig(fig, f"{plot_fname}.png", dpi=72)
        # vector
        genetools.plots.savefig(fig, f"{plot_fname}.pdf")

        plt.close(fig)

    # For each hue: plot, for points belonging to each hue value, the distribution of the distance between those points and their centroid
    for (hue_col, is_continuous) in colors_to_plot:
        if is_continuous:
            # Skip continuous hues
            continue

        # Find centroid and compute centroid->point distances in raw data space or in PCA space.
        for use_rep in ["X", "X_pca"]:
            # set up time course figure for this hue of interest. the figure is a single row of subplots (dim: # of timepoints)
            fig, axarr = plt.subplots(nrows=1, ncols=n_timepoints, figsize=(32, 6))

            # plot for each timepoint
            for timepoint_ix, ((timepoint, adata_transformed), ax) in enumerate(
                zip(timepoint_adatas_pcaed, axarr)
            ):
                plots.plot_distance_from_centroid_for_each_hue_value(
                    adata_transformed,
                    hue_col,
                    title=f"{timepoint} (reference)"
                    if timepoint == reference_timepoint
                    else timepoint,
                    ax=ax,
                    palette=palette_dict,
                    use_rep=use_rep,
                )

            fig.tight_layout()

            # save time course figure
            # convert plot name to valid filename
            plot_fname = f"{config.paths.output_dir}/{slugify(plot_name)}.distance_from_centroid.{slugify(hue_col)}.{use_rep}"
            # rasterized
            genetools.plots.savefig(fig, f"{plot_fname}.png", dpi=72)
            # vector
            genetools.plots.savefig(fig, f"{plot_fname}.pdf")

            plt.close(fig)
