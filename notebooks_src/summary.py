# -*- coding: utf-8 -*-
# %%

# %%
from summarynb import show
from IPython.display import display, Markdown
from covid_serology import config
from slugify import slugify

# %%
def summary(plot_name, hue_name, full=True):
    display(Markdown(f"## Scatterplot colored by {hue_name}"))
    show(
        f"{config.paths.output_dir}/{slugify(plot_name)}.pca.{slugify(hue_name)}.png",
        max_width=2000,
    )

    if full:
        # Don't stop early after just plotting the main hue_name's scatterplot
        display(Markdown(f"### Distance from centroid of each {hue_name} group"))
        for (use_rep, description) in zip(
            ["X", "X_pca"], ["raw data", "PCA (relative to reference timepoint)"]
        ):
            display(Markdown(f"#### In {description} space"))
            show(
                f"{config.paths.output_dir}/{slugify(plot_name)}.distance_from_centroid.{slugify(hue_name)}.{use_rep}.png",
                max_width=2000,
            )


# %%

# %% [markdown]
# # PCA colored by Exposure
#
# To characterize how spread out each group is, we also find each group’s centroid and plot the distribution of distances from all points in that group to their centroid. We measure the distance in raw data space and in PC1-PC2 space. Note that the PC space is created relative to a reference timepoint, while distances in raw data space in any timepoint are independent of other timepoints.
#

# %%

# %% [markdown]
# # Fig2C-D: CoV2 antigens (IgM, IgG, IgA) except N antigen

# %%
# Include Exposure Type in summary so you can see what the groups are
summary("cov2-all-except-n", "Exposure Type", full=False)

summary("cov2-all-except-n", "Exposure")

# %%

# %%

# %% [markdown]
# ---
#
# # Fig S5B: Time course PCA of variant ratios for `infection_cohort2` patients and Pfizer vaccinees AND variant infections
# NOTE: this variant infections PCA plot does not include breakthrough infections (i.e. post-vax)
# at each timepoint, only choose subjects who don't violate any p11 cutoffs at that timepoint.

# %%
# Include Exposure Type in summary so you can see what the groups are
summary("variants-with-variant-infections", "Exposure Type", full=False)

summary("variants-with-variant-infections", "Exposure")

# %% [markdown]
# ---
#
# # Day 90 variant ratios for `infection_cohort2` patients, Pfizer vaccinees, Mongolia vaccinees. Without vaccinees who are CoV2+

# %%
show(
    f"{config.paths.output_dir}/variant-ratios-all-groups-day90-week-7-and-later-3-months-exposure-type.pca.png",
    max_width=2000,
)

# %% [markdown]
# ---
#
# # Fig 5B: Imprinting index
# Below are the preference/imprinting levels for:
#
# - Delta infected
# - Wuhan vaccinated
# - Wuhan vaccinated then Delta infected

# %%
show(f"{config.paths.output_dir}/imprinting_level.log_rescaled_ratios.png")

# %% [markdown]
# #### Preliminary: Plot raw data, Wuhan vs variant
# The $y=x$ dotted line indicates Wuhan/variant ratios of 1. Above the line means ratio < 1, and below the line means ratio > 1.

# %%
# show(
#     f"{config.paths.output_dir}/imprinting_level.raw_values.png",
#     headers=["Raw values"],
#     max_width=2000,
# )
show(
    [
        f"{config.paths.output_dir}/imprinting_level.raw_values.{slugify(status)}.png"
        for status in [
            "Infected with Delta",
            "Vaccinated against Wuhan + Infected with Delta",
            "Vaccinated against Wuhan",
        ]
    ],
    headers=["Raw values"],
)

# %%

# %% [markdown]
# #### Plot ratios and transformed ratios
# The dotted lines indicate where Wuhan/variant preference is even (i.e. raw ratio = 1 or log ratio = 0)

# %%
show(
    [
        f"{config.paths.output_dir}/imprinting_level.raw_ratios.png",
        f"{config.paths.output_dir}/imprinting_level.log_ratios.png",
        f"{config.paths.output_dir}/imprinting_level.log_rescaled_ratios.png",
        f"{config.paths.output_dir}/imprinting_level.log_rescaled_ratios.densities.stacked.png",
    ],
    headers=[
        "Raw ratios",
        "Log transformed ratios",
        "Log transformed and rescaled ratios (imprinting amounts)",
        "Distribution of log transformed and rescaled ratios (imprinting amounts)",
    ],
)

# %%
