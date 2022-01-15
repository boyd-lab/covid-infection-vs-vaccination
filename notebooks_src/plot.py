# %% [markdown]
# # PCA timepoints independently, and PCA timepoints against a common reference timepoint

# %%

# %%
import numpy as np
import matplotlib.pyplot as plt

# %matplotlib inline
import pandas as pd
import seaborn as sns
import scanpy as sc
import genetools
from genetools.palette import HueValueStyle
from slugify import slugify
from covid_serology import config
from covid_serology.main import (
    transform_and_plot_pcas_multiple_timepoints_against_common_reference,
)
from covid_serology.pca import (
    pca_one_timepoint_against_reference,
    cleanup_adata,
)
from covid_serology.plots import plot_pca_one_hue_one_timepoint
from covid_serology.helpers import palette_dict

# %%
# Set fonts
from matplotlib import rc

rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"], "size": 12})


# %%
# # view the params with:
# sns.plotting_context("paper",font_scale=1.25)

sns.set_context("paper", font_scale=1.25)

# %%
sc.settings.verbosity = 3

# %%

# %%

# %% [markdown]
# ## Load datasets
# Timepoints are concatenated horizontally. Use "timepoint" annotation in adata.var to separate out.
#
# Rows are subjects. Columns are measurements defined as `(virus, target, variant_plate_type, antibody, timepoint)` tuples. Given a measurement at a timepoint: the same measurement at a different timepoint will be stored in a different column.

# %%

# %%

# %%
adata_all = sc.read(f"{config.paths.generated_data_dir}/combined_data.h5")
adata_all

# %%
adata_all.obs

# %%
adata_all.var

# %%


# %%

# %%
adata_variant_ratios = sc.read(f"{config.paths.generated_data_dir}/variant_ratios.h5")
adata_variant_ratios

# %%
adata_variant_ratios.obs["source_cohort"].value_counts()

# %%
# Remove any subjects that now have NaNs for all variables
# These will be any samples that have no variant ratio measurements
sc.pp.filter_cells(adata_variant_ratios, min_genes=1)

adata_variant_ratios.shape

# %%
adata_variant_ratios.obs["source_cohort"].value_counts()

# %%

# %%

# %%
adata_coronavirus_plate_subset = sc.read(
    f"{config.paths.generated_data_dir}/coronavirus_plate_only.subset.h5"
)
adata_coronavirus_plate_subset

# %%
adata_coronavirus_plate_subset.obs

# %%
adata_coronavirus_plate_subset.var

# %%
# Filter to Wuhan/non-variant only
adata_coronavirus_plate_subset = adata_coronavirus_plate_subset[
    :, adata_coronavirus_plate_subset.var["variant_plate_type"] == "Wuhan"
].copy()
adata_coronavirus_plate_subset.var

# %%
adata_coronavirus_plate_subset.obs["source_cohort"].value_counts()

# %%
adata_coronavirus_plate_subset.shape

# %%
# Remove any subjects that now have NaNs for all variables
# These will be any samples that have no standard coronavirus plate measurements, only variant ratios
sc.pp.filter_cells(adata_coronavirus_plate_subset, min_genes=1)

adata_coronavirus_plate_subset.shape

# %%
adata_coronavirus_plate_subset.obs["source_cohort"].value_counts()

# %%
# %%

# %%

# %% [markdown]
# ### Choose timepoints

# %%


# %%
all_timepoints = set(adata_coronavirus_plate_subset.var["timepoint"]).union(
    adata_variant_ratios.var["timepoint"]
)
all_timepoints
# %%
reference_timepoint_wuhan_comparison = config.reference_timepoint_wuhan_comparison
reference_timepoint_wuhan_comparison

# %%
reference_timepoint_variant_comparison = config.reference_timepoint_variant_comparison
reference_timepoint_variant_comparison

# %%
timepoint_sort_order = [
    "day 0 / pre-pandemic",
    "day 7 / week 1",
    "day 21 / weeks 2&3",
    "day 28 / week 4",
    "day 42 / weeks 5&6",
    "week 7 and later / 3 months",
    "day 210 / 7 months",
    "boostD1/2",
    "boostD21",
    "boostD7/8",
]
assert all(
    [timepoint in timepoint_sort_order for timepoint in all_timepoints]
), "Not all timepoints are represented in timepoint sort order"
timepoint_sort_order

# %%

# %%

# %% [markdown]
# ### Choose colors

# %%
discrete_colors_to_plot = [
    "Exposure",
    "Exposure Type",
]
# could add continuous hues here too

# (column name, continuous or not)
colors_to_plot = [(s, False) for s in discrete_colors_to_plot]
colors_to_plot

# %%
palette_dict.keys()

# %%

# %%
# visualize chosen colors
# first cast any HueValueStyle objects in the dict to be color strings only
palette_dict_colors_only = HueValueStyle.huestyles_to_colors_dict(palette_dict)

# modified palplot to have labels: https://stackoverflow.com/a/64492813 and https://github.com/mwaskom/seaborn/blob/cef0a2d5e86477a80659898e66bc2295886bf917/seaborn/miscplot.py#L9

import matplotlib.ticker as ticker
import matplotlib as mpl

n = len(palette_dict_colors_only)
size = 1

f, ax = plt.subplots(1, 1, figsize=(n * size, size))
ax.imshow(
    np.arange(n).reshape(1, n),
    cmap=mpl.colors.ListedColormap(list(palette_dict_colors_only.values())),
    interpolation="nearest",
    aspect="auto",
)
ax.set_xticks(np.arange(n) - 0.5)
ax.set_yticks([-0.5, 0.5])
# Ensure nice border between colors
ax.set_xticklabels(["" for _ in range(n)])

# Add labels to x-axis:
ax.set_xticks(np.arange(n))
ax.set_xticklabels(palette_dict_colors_only.keys(), ha="center")

# The proper way to set no ticks on y-axis
ax.yaxis.set_major_locator(ticker.NullLocator())

genetools.plots.wrap_tick_labels(ax, wrap_amount=10)


# %%

# %% [markdown]
# ## CoV2 antigens, IgG, IgM, and IgA only, except N antigen
#
# Using subset of data:
#
# - `infection_cohort1`
# - Pfizer vaccine data, original version (new version no longer has IgM or IgA data, only IgG)

# %%
adata_coronavirus_plate_subset

# %%
adata_coronavirus_plate_subset.obs["source_cohort"].value_counts()

# %%
adata_coronavirus_plate_subset.var

# %%
_ = transform_and_plot_pcas_multiple_timepoints_against_common_reference(
    adata=adata_coronavirus_plate_subset[
        :,
        (adata_coronavirus_plate_subset.var["virus"] == "Wuhan")
        & (adata_coronavirus_plate_subset.var["target"] != "N"),
    ],
    plot_name="cov2_all_except_N",
    palette_dict=palette_dict,
    colors_to_plot=colors_to_plot,
    reference_timepoint=reference_timepoint_wuhan_comparison,
    timepoint_sort_order=timepoint_sort_order,
)

# %%

# %%

# %% [markdown]
# ## Time course PCA of variant ratios for `infection_cohort2` patients and Pfizer vaccinees AND variant infections
# at each timepoint, only choose subjects who don't violate any p11 cutoffs at that timepoint.


# %%
adata_variant_ratios.obs["source_cohort"].value_counts()

# %%

# %%
# subset to cohorts
adata_variant_ratios_time_course_with_variant_infections = adata_variant_ratios[
    adata_variant_ratios.obs["source_cohort"].isin(
        ["infection_cohort2_all", "Pfizer_all", "variant_infections"]
    )
].copy()
adata_variant_ratios_time_course_with_variant_infections.obs["Status"].value_counts()

# %%
adata_variant_ratios_time_course_with_variant_infections

# %%

# %%
# also subset to not being a breakthrough infection
adata_variant_ratios_time_course_with_variant_infections = (
    adata_variant_ratios_time_course_with_variant_infections[
        adata_variant_ratios_time_course_with_variant_infections.obs["Status"]
        != "Variant Infection - Delta - mRNA vaccinated"
    ].copy()
)
adata_variant_ratios_time_course_with_variant_infections.obs["Status"].value_counts()

# %%
adata_variant_ratios_time_course_with_variant_infections

# %%

# %%
adata_variant_ratios_time_course_with_variant_infections.obs[
    "any_p11_cutoffs_violated"
].value_counts()

# %%
# timepoint specific cutoffs violated columns:
adata_variant_ratios_time_course_with_variant_infections.obs[
    "any_p11_cutoffs_violated:day 28 / week 4"
].value_counts()

# %%
adata_variant_ratios_time_course_with_variant_infections.var

# %%

# %%
adata_variant_ratios_time_course_with_variant_infections_figure_data = (
    transform_and_plot_pcas_multiple_timepoints_against_common_reference(
        adata=adata_variant_ratios_time_course_with_variant_infections,
        plot_name="variants_with_variant_infections",
        palette_dict=palette_dict,
        colors_to_plot=colors_to_plot,
        reference_timepoint=reference_timepoint_variant_comparison,
        enforce_cutoffs=True,
        timepoint_sort_order=timepoint_sort_order,
    )
)

# %%

# %% [markdown]
# ### Export exactly the samples that are in the figure (i.e. passed cutoffs)

# %%
adata_variant_ratios_time_course_with_variant_infections_figure_data

# %%
# reformat as dict
adata_variant_ratios_time_course_with_variant_infections_figure_data = {
    timepoint: timepoint_adata_pcaed
    for (
        timepoint,
        timepoint_adata_pcaed,
    ) in adata_variant_ratios_time_course_with_variant_infections_figure_data
}
adata_variant_ratios_time_course_with_variant_infections_figure_data.keys()

# %%
# get raw data (pre-transforms) for reference timepoint
fig_s5b_adata = adata_variant_ratios_time_course_with_variant_infections_figure_data[
    reference_timepoint_variant_comparison
].raw.to_adata()
fig_s5b_adata

# %%

# %% [markdown]
# Alternatively you could use:
#
# ```
# fig_s5b_adata = adata_variant_ratios_time_course_with_variant_infections[
#     adata_variant_ratios_time_course_with_variant_infections.obs[
#         f"any_p11_cutoffs_violated:{reference_timepoint_variant_comparison}"
#     ]
#     == False,
#     adata_variant_ratios_time_course_with_variant_infections.var["timepoint"]
#     == reference_timepoint_variant_comparison,
# ]
# fig_s5b_adata = cleanup_adata(fig_s5b_adata)
# ```
#
# But that's more error-prone. Let's let the main plotting code handle that logic for us
#
#
# Continue with the export:

# %%

# %%
fig_s5b_df = fig_s5b_adata.to_df()
fig_s5b_df

# %%
fig_s5b_df.isna().any()

# %%
# Reshape
fig_s5b_df = pd.melt(
    fig_s5b_df.reset_index(),
    id_vars=["index"],
    var_name="var_name",
    value_name="ratio_value",
).rename(columns={"index": "participant"})
fig_s5b_df

# %%
fig_s5b_df["var_name"].nunique()

# %%
fig_s5b_adata.obs["Exposure"]

# %%
# Merge in Exposure
fig_s5b_df = pd.merge(
    fig_s5b_df,
    fig_s5b_adata.obs["Exposure"].rename("exposure"),
    left_on="participant",
    right_index=True,
    validate="m:1",
    how="left",
)
assert not fig_s5b_df["exposure"].isna().any(), "merge failed"
fig_s5b_df

# %%
fig_s5b_adata.var

# %%
# Merge in var info
fig_s5b_df = pd.merge(
    fig_s5b_df,
    fig_s5b_adata.var.drop("n_cells", axis=1),
    left_on="var_name",
    right_index=True,
    validate="m:1",
    how="left",
)
assert not fig_s5b_df["numerator"].isna().any(), "merge failed"
fig_s5b_df

# %%
# Reorder and simplify
fig_s5b_df = fig_s5b_df.rename(columns={"virus": "ratio_denominator_strain"}).assign(
    ratio_numerator_strain="Wuhan"
)
fig_s5b_df = fig_s5b_df[
    [
        "participant",
        "exposure",
        "timepoint",
        "ratio_numerator_strain",
        "ratio_denominator_strain",
        "ratio_value",
    ]
]
fig_s5b_df = fig_s5b_df.sort_values(["participant", "ratio_denominator_strain"])
fig_s5b_df

# %%
fig_s5b_df.to_csv(
    f"{config.paths.output_dir}/fig_s5b_data_export.tsv", sep="\t", index=None
)

# %%

# %%

# %%

# %% [markdown]
# ## Day 90 variant ratios for `infection_cohort2` patients, Pfizer vaccinees, Mongolia vaccinees. Without vaccinees who are CoV2+

# %%
adata_variant_ratios.obs["source_cohort"].value_counts()

# %%

# %%
# subset cohorts
adata_variant_ratios_day90 = adata_variant_ratios[
    adata_variant_ratios.obs["source_cohort"] != "infection_cohort1"
].copy()
adata_variant_ratios_day90.obs["Exposure Type"].value_counts()

# %%
adata_variant_ratios_day90

# %%

# %%
# filter out breakthrough infections or past infections before vaccination

adata_variant_ratios_day90 = adata_variant_ratios_day90[
    ~adata_variant_ratios_day90.obs["Exposure Type"].isin(
        [
            "Pfizer-Pfizer (Stanford), CoV2+",
            "Variant Infection - Delta - mRNA vaccinated",
        ]
    )
].copy()

adata_variant_ratios_day90.obs["Exposure Type"].value_counts()

# %%
adata_variant_ratios_day90

# %%

# %%
adata_variant_ratios_day90.var

# %%
adata_variant_ratios_day90.var["timepoint"].value_counts()

# %%
timepoint = "week 7 and later / 3 months"
plot_name = "variant_ratios_all_groups_day90"

# convert plot name to valid filename
output_prefix = slugify(f"{plot_name}.{timepoint}")

adata_subset = adata_variant_ratios_day90[
    adata_variant_ratios_day90.obs[f"any_p11_cutoffs_violated:{timepoint}"] == False,
    adata_variant_ratios_day90.var["timepoint"] == timepoint,
]
adata_subset = cleanup_adata(adata_subset)
adata_transformed = pca_one_timepoint_against_reference(
    adata_subset, output_prefix=output_prefix, transformer=None
)


# plot PCA
for (hue_col, is_continuous) in colors_to_plot:
    fig, ax = plt.subplots(figsize=(6, 6))
    plot_pca_one_hue_one_timepoint(
        adata=adata_transformed,
        ax=ax,
        hue_key=hue_col,
        title=f"{hue_col} at {timepoint}",
        enable_legends=True,
        palette=palette_dict,
        continuous_hue=is_continuous,
    )
    plot_fname_hue = slugify(f"{plot_name}.{timepoint}.{hue_col}")
    # save rasterized
    genetools.plots.savefig(
        fig, f"{config.paths.output_dir}/{plot_fname_hue}.pca.png", dpi=72
    )
    # save vector
    genetools.plots.savefig(fig, f"{config.paths.output_dir}/{plot_fname_hue}.pca.pdf")

# %%

# %%

# %%

# %% [markdown]
# ## Imprinting amount

# %% [markdown]
# ### Get data

# %%
adata_variant_ratios.obs["source_cohort"].value_counts()

# %%
adata_imprinting = adata_variant_ratios[
    adata_variant_ratios.obs["source_cohort"].isin(["Pfizer_all", "variant_infections"])
].copy()
adata_imprinting.obs["Status"].value_counts()

# %%
adata_imprinting = adata_imprinting[
    adata_imprinting.obs["Status"] != "Variant Infection - Alpha"
].copy()
adata_imprinting.obs["Status"].value_counts()

# %%
adata_imprinting.obs["Exposure Type"].value_counts()

# %%
adata_imprinting = adata_imprinting[
    adata_imprinting.obs["Exposure Type"] != "Pfizer-Pfizer (Stanford), CoV2+"
].copy()
adata_imprinting.obs["Status"].value_counts()

# %%
adata_imprinting.obs["Exposure Type"].value_counts()

# %%
adata_imprinting

# %%

# %%
adata_imprinting = adata_imprinting[
    :,
    (adata_imprinting.var["virus"] == "Delta")
    & (adata_imprinting.var["timepoint"] == reference_timepoint_variant_comparison),
]
adata_imprinting

# %%
adata_imprinting = adata_imprinting[
    ~adata_imprinting.obs[
        f"any_p11_cutoffs_violated:{reference_timepoint_variant_comparison}"
    ]
]

# %%
assert adata_imprinting.var.shape[0] == 1
adata_imprinting.var

# %%

# %%

# %%

# %%
imprinting_df = adata_imprinting.obs[["Status"]].assign(ratio=adata_imprinting.X)
imprinting_df

# %%
imprinting_df.dropna(inplace=True)
imprinting_df

# %%
imprinting_df["Status"].value_counts()

# %%
imprinting_df["Status"] = imprinting_df["Status"].replace(
    {
        "Vaccinee": "Vaccinated against Wuhan",
        "Variant Infection - Delta": "Infected with Delta",
        "Variant Infection - Delta - mRNA vaccinated": "Vaccinated against Wuhan + Infected with Delta",
    }
)

# %%
imprinting_df["Status"].value_counts()

# %%

# %% [markdown]
# ### Plot raw data

# %%
imprinting_df

# %%
# get raw values
numerator_varname = adata_imprinting.var.squeeze()["numerator"]
denominator_varname = adata_imprinting.var.squeeze()["denominator"]

imprinting_df = genetools.helpers.merge_into_left(
    imprinting_df, adata_all[:, [numerator_varname, denominator_varname]].to_df()
)
assert (
    not imprinting_df[[numerator_varname, denominator_varname]].isna().any().any()
), "some raw values missing"
imprinting_df

# %%
g = sns.FacetGrid(
    imprinting_df,
    col="Status",
    hue="Status",
    height=4,
    sharey=False,
    sharex=False,
)
g.map(sns.scatterplot, numerator_varname, denominator_varname, alpha=0.7)
for ax in g.axes_dict.values():
    ax.axline((0, 0), slope=1, color=".2", linestyle="--", zorder=0)

for ax in g.axes.flat:
    # Make square: https://stackoverflow.com/q/25497402/130164
    # Get min and max of both axes
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),
        np.max([ax.get_xlim(), ax.get_ylim()]),
    ]
    ax.set_aspect("equal")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
g.add_legend()

genetools.plots.savefig(
    g.fig, f"{config.paths.output_dir}/imprinting_level.raw_values.png", dpi=72
)

# %%

# %%
with sns.axes_style("white"):
    for (status, grp), color in zip(
        imprinting_df.groupby("Status", sort=False), sns.color_palette()
    ):
        g = sns.jointplot(
            data=grp,
            x=numerator_varname,
            y=denominator_varname,
            color=color,
        )
        ax = g.ax_joint
        ax.axline((0, 0), slope=1, color=".2", linestyle="--", zorder=0)

        # Make square: https://stackoverflow.com/q/25497402/130164
        # Get min and max of both axes
        lims = [
            np.min([ax.get_xlim(), ax.get_ylim()]),
            np.max([ax.get_xlim(), ax.get_ylim()]),
        ]
        ax.set_aspect("equal")
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        g.fig.subplots_adjust(top=0.9)
        g.fig.suptitle(status)
        genetools.plots.savefig(
            g.fig,
            f"{config.paths.output_dir}/imprinting_level.raw_values.{slugify(status)}.png",
            dpi=72,
        )

# %%

# %% [markdown]
# ### Plot raw ratios

# %%
fig, ax = plt.subplots()
sns.boxplot(data=imprinting_df, x="Status", y="ratio", ax=ax)
sns.swarmplot(data=imprinting_df, x="Status", y="ratio", color=".25", alpha=0.7, ax=ax)

# Even preferences line
plt.axhline(y=1, color="k", zorder=-5, linestyle="dashed")

plt.ylabel("Wuhan / Delta ratio")

# wrap tick labels
genetools.plots.wrap_tick_labels(ax, wrap_x_axis=True, wrap_y_axis=False)

genetools.plots.savefig(
    fig, f"{config.paths.output_dir}/imprinting_level.raw_ratios.png", dpi=72
)

# %%

# %% [markdown]
# ### Plot log transformed ratio

# %%
imprinting_df["ratio_log"] = np.log(imprinting_df["ratio"])
imprinting_df

# %%
fig, ax = plt.subplots()

sns.boxplot(data=imprinting_df, x="Status", y="ratio_log", ax=ax)
sns.swarmplot(
    data=imprinting_df, x="Status", y="ratio_log", color=".25", alpha=0.7, ax=ax
)

# even preference line
plt.axhline(y=0, color="k", zorder=-5, linestyle="dashed")

# wrap tick labels
genetools.plots.wrap_tick_labels(ax, wrap_x_axis=True, wrap_y_axis=False)

plt.ylabel("Wuhan / Delta ratio, log transformed")

genetools.plots.savefig(
    fig, f"{config.paths.output_dir}/imprinting_level.log_ratios.png", dpi=72
)

# %%
imprinting_df["ratio_log"].describe()

# %%

# %% [markdown]
# ### Rescale log-transformed ratio, piecewise: -100%-0 and 0-100%

# %%

# %%
# ratio_log_rescaled
imprinting_df["Imprinting level"] = imprinting_df["ratio_log"].copy()
maximum_ratio_log = np.abs(imprinting_df["ratio_log"].max())
minimum_ratio_log = np.abs(imprinting_df["ratio_log"].min())
imprinting_df.loc[imprinting_df["ratio_log"] > 0, "Imprinting level"] = (
    imprinting_df.loc[imprinting_df["ratio_log"] > 0, "ratio_log"] / maximum_ratio_log
)
imprinting_df.loc[imprinting_df["ratio_log"] < 0, "Imprinting level"] = (
    imprinting_df.loc[imprinting_df["ratio_log"] < 0, "ratio_log"] / minimum_ratio_log
)
imprinting_df["Imprinting level"].describe()

# %%

# %%
fig, ax = plt.subplots()
sns.boxplot(data=imprinting_df, x="Status", y="Imprinting level", ax=ax)
sns.swarmplot(
    data=imprinting_df, x="Status", y="Imprinting level", color=".25", alpha=0.7, ax=ax
)
plt.axhline(y=0, color="k", zorder=-5, linestyle="dashed")
plt.yticks(
    ticks=[-1, -0.5, 0, 0.5, 1],
    labels=[
        "100% Delta preference",
        "50%",
        "Even preference",
        "50%",
        "100% Wuhan preference",
    ],
)

# Add sample size to labels
ax.set_xticklabels(
    genetools.plots.add_sample_size_to_labels(
        labels=ax.get_xticklabels(), data=imprinting_df, hue_key="Status"
    )
)

# wrap tick labels
genetools.plots.wrap_tick_labels(ax, wrap_x_axis=True, wrap_y_axis=True)


# save rasterized
genetools.plots.savefig(
    fig, f"{config.paths.output_dir}/imprinting_level.log_rescaled_ratios.png", dpi=72
)
# save vector
genetools.plots.savefig(
    fig, f"{config.paths.output_dir}/imprinting_level.log_rescaled_ratios.pdf"
)

# %%

# %%
# Same but densities, overlapping

fig, ax = plt.subplots(figsize=(8, 4))
ax = sns.kdeplot(
    data=imprinting_df, x="Imprinting level", hue="Status", ax=ax, legend=True
)
# Move legend -- need special way to extract legend for histplot or kdeplot: https://github.com/mwaskom/seaborn/issues/2280#issuecomment-899101193
# place legend outside figure
sns.move_legend(
    ax,
    bbox_to_anchor=(1.05, 0.5),
    loc="center left",
    borderaxespad=0.0,
    # no border
    frameon=False,
    # transparent background
    framealpha=0.0,
)


plt.axvline(x=0, color="k", zorder=-5, linestyle="dashed")
plt.xticks(
    ticks=[-1, -0.5, 0, 0.5, 1],
    labels=[
        "100% Delta preference",
        "50%",
        "Even preference",
        "50%",
        "100% Wuhan preference",
    ],
)

# wrap tick labels
genetools.plots.wrap_tick_labels(ax, wrap_x_axis=True, wrap_y_axis=False)

genetools.plots.savefig(
    fig,
    f"{config.paths.output_dir}/imprinting_level.log_rescaled_ratios.densities.overlap.png",
    dpi=72,
)

# %%
# Same but densities, stacked

g = sns.FacetGrid(
    imprinting_df, row="Status", hue="Status", sharey=False, sharex=True, aspect=2
)
g.map(sns.kdeplot, "Imprinting level")

for ax in g.axes.flat:
    ax.axvline(x=0, color=".2", zorder=-5, linestyle="dashed", alpha=0.7)

    ax.set_xticks(ticks=[-1, 0, 1], minor=False)
    ax.set_xticklabels(
        [
            "100% Delta preference",
            "Even preference",
            "100% Wuhan preference",
        ],
        minor=False,
    )

    ax.set_xticks(ticks=[-0.5, 0.5], minor=True)
    ax.set_xticklabels(
        [
            "50%",
            "50%",
        ],
        minor=True,
    )

    # wrap tick labels
    genetools.plots.wrap_tick_labels(ax, wrap_x_axis=True, wrap_y_axis=False)


# save rasterized
genetools.plots.savefig(
    g.fig,
    f"{config.paths.output_dir}/imprinting_level.log_rescaled_ratios.densities.stacked.png",
    dpi=72,
)
# save vector
genetools.plots.savefig(
    g.fig,
    f"{config.paths.output_dir}/imprinting_level.log_rescaled_ratios.densities.stacked.pdf",
)

# %%

# %%

# %%

# %%
