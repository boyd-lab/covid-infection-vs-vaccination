# %%

# %%
import numpy as np
import matplotlib.pyplot as plt

# %matplotlib inline
import pandas as pd
import seaborn as sns
import anndata
import scanpy as sc
import genetools
from covid_serology import config
from numpy.testing import assert_array_equal

# %%

# %% [markdown]
# # Load data

# %%
# patients have multiple timepoints
infection_df = pd.read_csv(f"{config.paths.data_dir}/infection_cohort2.coronavirus.csv")
infection_df

# %%

# %%
infection_participants = pd.read_csv(
    f"{config.paths.data_dir}/infection_cohort2.demographics.csv"
)
infection_participants

# %%
# one row per patient in demographics table
assert all(infection_participants["Record.ID"].value_counts() == 1)

# %%
infection_participants.columns

# %%

# %% [markdown]
# # infection data reshape

# %%
infection_df.columns

# %%
# only one sample per patient per timepoint
assert all(infection_df.groupby(["sample_id", "timepoint"]).size() == 1)

# %%
# how many total samples per patient? many have only a single sample:
infection_df["sample_id"].value_counts().value_counts()

# %%

# %%
infection_df = pd.melt(
    infection_df,
    id_vars=["sample_id", "timepoint", "days_post_pcr"],
    var_name="measurement",
)
infection_df

# %%
# just in case, convert to float and switch non-numeric values to nan
infection_df["value"] = pd.to_numeric(infection_df["value"], errors="coerce")
infection_df.dtypes

# %%
infection_df["value"].isna().value_counts()

# %%
infection_df[infection_df["value"].isna()]

# %%

# %%
infection_df["measurement"].unique()

# %%
infection_df = infection_df[infection_df["measurement"].str.contains("AU")].copy()
infection_df["measurement"].value_counts()

# %%
infection_df["measurement"] = infection_df["measurement"].str.replace("_AU", "")
infection_df["measurement"].value_counts()

# %%

# %%
# extract parts of measurement column
# all coronavirus plate measurements are IgG only
measurement_parts = (
    infection_df["measurement"]
    .str.split("_", expand=True)
    .rename(columns={0: "virus", 1: "target"})
    .assign(variant_plate_type="Wuhan", antibody="IgG")
    .apply(lambda col: col.str.strip())
)
measurement_parts

# %%
measurement_parts["virus"].value_counts()

# %%
measurement_parts["target"].value_counts()

# %%
measurement_parts["antibody"].value_counts()

# %%
measurement_parts["variant_plate_type"].value_counts()

# %%

# %%
infection_df = pd.concat([infection_df, measurement_parts], axis=1)
infection_df

# %%

# %%
infection_df["days_post_pcr"].hist()

# %%
infection_df["timepoint"].value_counts()

# %%
# timepoint label map
map_infection_to_global_timepoint_labels = {
    1: "day 21 / weeks 2&3",
    2: "day 28 / week 4",
    3: "week 7 and later / 3 months",
    4: "day 210 / 7 months",
}
assert all(
    tp in map_infection_to_global_timepoint_labels.keys()
    for tp in infection_df["timepoint"].unique()
)

# %%
infection_df


# %%

# %%

# %%
def process_infection_timepoint(df_partial, timepoint):
    # at a given time point: only one measurement per patient-virus-target combo
    assert all(
        df_partial.groupby(
            ["sample_id", "virus", "target", "variant_plate_type", "antibody"]
        ).size()
        == 1
    )

    # unmelt into matrix
    infection_df_pivot = pd.pivot(
        df_partial,
        index="sample_id",
        columns=[
            "virus",
            "target",
            "variant_plate_type",
            "antibody",
        ],
        values="value",
    )

    ## set column names
    variable_info = infection_df_pivot.columns.to_frame().reset_index(drop=True)
    # create combined name
    variable_info["timepoint"] = timepoint
    variable_info["combined_name"] = variable_info.apply("_".join, axis=1)
    variable_info = variable_info.set_index("combined_name")

    # set var names
    infection_df_pivot.columns = variable_info.index.copy()

    # drop patients with any NaNs in this timepoint
    infection_df_pivot = infection_df_pivot.dropna(how="any")
    assert not infection_df_pivot.isna().any().any()

    return infection_df_pivot, variable_info


# %%

# %%
X_partial = []
var_partial = []

for infection_timepoint in infection_df["timepoint"].unique():
    associated_global_timepoint_label = map_infection_to_global_timepoint_labels[
        infection_timepoint
    ]

    print(infection_timepoint, "->", associated_global_timepoint_label)

    df_partial = infection_df[infection_df["timepoint"] == infection_timepoint]
    infection_df_pivot, variable_info = process_infection_timepoint(
        df_partial, associated_global_timepoint_label
    )

    X_partial.append(infection_df_pivot)
    var_partial.append(variable_info)

infection_df_pivot = pd.concat(X_partial, axis=1)
variable_info = pd.concat(var_partial, axis=0)

# %%
# note: there are NaNs - patients don't have entries for all timepoints
infection_df_pivot

# %%
variable_info

# %%

# %%

# %%
infection_participants

# %%
# attach status
infection_participants = infection_participants.set_index("Record.ID")
# set index to str for anndata: anndata requires string indices
infection_participants.index = infection_participants.index.astype(str)

# %%
infection_participants["Status"].value_counts()

# %%
# reorder participants info to match
infection_participants = infection_participants.loc[infection_df_pivot.index]
infection_participants

# %%
infection_participants.columns

# %%

# %%

# %%
adata_infection = anndata.AnnData(
    X=infection_df_pivot, obs=infection_participants, var=variable_info
)
adata_infection

# %%
adata_infection.var

# %%
adata_infection.obs["Status"].value_counts()

# %%

# %%
adata_infection.write(
    f"{config.paths.generated_data_dir}/partial.infection_cohort2.coronavirus_plate.h5"
)

# %%

# %%
# confirm multiple timepoints are kept for a patients with multiple timepoint
adata_infection.to_df().loc[infection_df["sample_id"].value_counts().idxmax()]

# %%
