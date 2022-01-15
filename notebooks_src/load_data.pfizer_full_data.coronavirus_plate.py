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
vaccine_participants = pd.read_csv(
    f"{config.paths.data_dir}/pfizer_demographics.csv"
).dropna(how="all")
vaccine_participants

# %%
vaccine_df = pd.read_csv(f"{config.paths.data_dir}/pfizer_coronavirus.csv")
vaccine_df

# %%
# only one sample per patient per timepoint
assert all(vaccine_df.groupby(["participant", "timepoint"]).size() == 1)

# %%

# %%

# %%

# %% [markdown]
# # Reshape vaccine data

# %%
vaccine_df.columns

# %%

# %%
vaccine_df = pd.melt(
    vaccine_df, id_vars=["participant", "timepoint"], var_name="measurement"
)
vaccine_df

# %%

# %%
# just in case, convert to float and switch non-numeric values to nan
vaccine_df["value"] = pd.to_numeric(vaccine_df["value"], errors="coerce")
vaccine_df.dtypes

# %%
vaccine_df["value"].isna().value_counts()

# %%
vaccine_df[vaccine_df["value"].isna()]

# %%

# %%
vaccine_df["measurement"].unique()

# %%

# %%
vaccine_df["measurement"] = vaccine_df["measurement"].str.replace("_AU", "")
vaccine_df["measurement"].value_counts()

# %%
# extract parts of measurement column
# all coronavirus plate measurements are IgG only
measurement_parts = (
    vaccine_df["measurement"]
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
measurement_parts["variant_plate_type"].value_counts()

# %%

# %%

# %%
vaccine_df = pd.concat([vaccine_df, measurement_parts], axis=1)
vaccine_df

# %%
vaccine_df["timepoint"].value_counts()

# %%
# timepoint label map
map_vaccine_to_global_timepoint_labels = {
    "D0": "day 0 / pre-pandemic",
    "D7": "day 7 / week 1",
    "D21": "day 21 / weeks 2&3",
    "D28": "day 28 / week 4",
    "D42": "day 42 / weeks 5&6",
    "D90": "week 7 and later / 3 months",
    "D210": "day 210 / 7 months",
    "boostD1/2": "boostD1/2",
    "boostD7/8": "boostD7/8",
    "boostD21": "boostD21",
}
assert all(
    tp in map_vaccine_to_global_timepoint_labels.keys()
    for tp in vaccine_df["timepoint"].unique()
)


# %%

# %%
def process_vaccine_timepoint(df_partial, timepoint):
    # at a given time point: only one measurement per patient-virus-target combo
    assert all(
        df_partial.groupby(
            ["participant", "virus", "target", "variant_plate_type", "antibody"]
        ).size()
        == 1
    )

    # unmelt
    vaccine_df_pivot = pd.pivot(
        df_partial,
        index="participant",
        columns=[
            "virus",
            "target",
            "variant_plate_type",
            "antibody",
        ],
        values="value",
    )

    ## set column names
    variable_info = vaccine_df_pivot.columns.to_frame().reset_index(drop=True)
    # create combined name
    variable_info["timepoint"] = timepoint
    variable_info["combined_name"] = variable_info.apply("_".join, axis=1)
    variable_info = variable_info.set_index("combined_name")

    # set var names
    vaccine_df_pivot.columns = variable_info.index.copy()

    # drop patients with any NaNs in this timepoint
    vaccine_df_pivot = vaccine_df_pivot.dropna(how="any")
    assert not vaccine_df_pivot.isna().any().any()

    return vaccine_df_pivot, variable_info


# %%

# %%
X_partial = []
var_partial = []
for vaccine_timepoint in vaccine_df["timepoint"].unique():
    associated_global_timepoint_label = map_vaccine_to_global_timepoint_labels[
        vaccine_timepoint
    ]
    print(vaccine_timepoint, "->", associated_global_timepoint_label)
    df_partial = vaccine_df[vaccine_df["timepoint"] == vaccine_timepoint]
    vaccine_df_pivot, variable_info = process_vaccine_timepoint(
        df_partial, associated_global_timepoint_label
    )
    X_partial.append(vaccine_df_pivot)
    var_partial.append(variable_info)
vaccine_df_pivot = pd.concat(X_partial, axis=1)
variable_info = pd.concat(var_partial, axis=0)

# %%
# note: there are NaNs - patients don't have entries for all timepoints
vaccine_df_pivot

# %%
variable_info

# %%

# %%
# attach status

vaccine_participants = vaccine_participants.set_index("PID")
vaccine_participants["Status"] = "Vaccinee"
# reorder
vaccine_participants = vaccine_participants.loc[vaccine_df_pivot.index]

# confirm same order
assert_array_equal(vaccine_participants.index, vaccine_df_pivot.index)

vaccine_participants

# %%
# anndata requires string indices
vaccine_participants.index = vaccine_participants.index.astype(str)
vaccine_df_pivot.index = vaccine_df_pivot.index.astype(str)

# %%

# %%
adata_vaccine = anndata.AnnData(
    X=vaccine_df_pivot, obs=vaccine_participants, var=variable_info
)
adata_vaccine

# %%
adata_vaccine.var

# %%
adata_vaccine.obs["Status"].value_counts()

# %%

# %%
adata_vaccine.write(
    f"{config.paths.generated_data_dir}/partial.pfizer_vaccine.coronavirus_plate.h5"
)

# %%

# %%

# %%
