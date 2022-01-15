# %%
import numpy as np
import matplotlib.pyplot as plt

# %matplotlib inline
import pandas as pd
import seaborn as sns
import anndata
import scanpy as sc
import genetools
from covid_serology import config, helpers
from numpy.testing import assert_array_equal

# %%

# %% [markdown]
# # Load data

# %%
mongolia_df = pd.read_csv(f"{config.paths.data_dir}/mongolia_mastersheet.csv")
mongolia_df

# %%
# Exclude vaccinees who also tested positive
# TODO: keep as a separate group

# use nucleocapsid cutoff, for all but sinopharm (which has killed virus)
nucleocapsid_cutoff_for_previous_infection = 3743.2
mongolia_df_without_sinopharm = mongolia_df[
    ~(
        (mongolia_df["first_dose_vaccine"] == "Sinopharm")
        | (mongolia_df["second_dose_vaccine"] == "Sinopharm")
    )
]
assert mongolia_df_without_sinopharm.shape[0] < mongolia_df.shape[0]
mongolia_df["evidence_of_previous_infection"] = (
    mongolia_df_without_sinopharm["Nucleocapsid"] > 3743.2
)
mongolia_df["evidence_of_previous_infection"]

# %%
mongolia_df["evidence_of_previous_infection"].value_counts(), mongolia_df[
    "evidence_of_previous_infection"
].isna().sum()

# %%
# for sinopharm, fillna with whether clinical record of positive test date exists
mongolia_df["evidence_of_previous_infection"] = mongolia_df[
    "evidence_of_previous_infection"
].fillna(~mongolia_df["tested_pos"].isna())

# all should now have a yes or no answer for whether previously infected
assert not mongolia_df["evidence_of_previous_infection"].isna().any()

mongolia_df["evidence_of_previous_infection"]

# %%
mongolia_df["evidence_of_previous_infection"].value_counts()

# %%
# Exclude vaccinees who also tested positive
# TODO: keep as a separate group
mongolia_df = mongolia_df[~mongolia_df["evidence_of_previous_infection"]]
mongolia_df

# %%

# %%
mongolia_df.columns

# %%
measurement_cols = mongolia_df.columns[mongolia_df.columns.str.contains("AU")]
helpers.confirm_all_measurement_columns_are_present(measurement_cols)
measurement_cols

# %%

# %%
mongolia_df["collection_after_vacc"].hist()

# %%
# no need to take means, already taken
assert all(mongolia_df.groupby(["sample_id", "collection_after_vacc"]).size()) == 1

# %%
# only one sample per patient
assert all(mongolia_df["sample_id"].value_counts() == 1)

# %%

# %% [markdown]
# # reformat

# %%

# %%
mongolia_df_obs = mongolia_df[
    [
        "sample_id",
        "collection_after_vacc",
        "first_dose_vaccine",
        "second_dose_vaccine",
    ]
].copy()

mongolia_df_obs.rename(columns={"collection_after_vacc": "Timepoint"}, inplace=True)

mongolia_df_obs["Status"] = (
    mongolia_df_obs["first_dose_vaccine"] + "-" + mongolia_df_obs["second_dose_vaccine"]
)

# capitalize columns for consistency - only the first letter
# and set index
mongolia_df_obs = mongolia_df_obs.rename(
    columns=lambda s: s[0].upper() + s[1:]
).set_index("Sample_id")

# anndata wants string index
mongolia_df_obs.index = mongolia_df_obs.index.astype(str)

mongolia_df_obs

# %%
# extract AU measurement cols, and set index to match obs
mongolia_df_X = (
    mongolia_df[measurement_cols]
    .rename(columns=lambda col: col.replace("_AU", ""))
    .set_index(mongolia_df_obs.index)
)
mongolia_df_X

# %%

# %%

# %%
adata = anndata.AnnData(X=mongolia_df_X, obs=mongolia_df_obs)
adata

# %%
adata.obs

# %%
adata.var

# %%
adata.var["virus"] = adata.var_names
adata.var["target"] = "RBD"
adata.var["variant_plate_type"] = "Variant"
adata.var["antibody"] = "IgG"
adata.var

# %%
adata.X

# %%
adata

# %%

# %%
adata.obs["Timepoint"].hist()

# %%
adata = adata[(adata.obs["Timepoint"] >= 60) & (adata.obs["Timepoint"] < 120)].copy()
adata

# %%
adata.var["timepoint"] = "week 7 and later / 3 months"
adata.var

# %%
# create combined name
adata.var["combined_name"] = adata.var.apply("_".join, axis=1)
adata.var = adata.var.set_index("combined_name")

# %%
adata.var

# %%
adata.obs["Status"].value_counts()

# %%

# %%

# %%

# %%
adata.write(f"{config.paths.generated_data_dir}/partial.mongolia_vaccines.h5")

# %%

# %%

# %%
