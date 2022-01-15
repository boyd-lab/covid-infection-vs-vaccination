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

# %%

# %% [markdown]
# # Load data

# %%
infection_df = pd.read_csv(
    f"{config.paths.data_dir}/variant_infections_seropositives.plate11.csv"
)
infection_df.head()

# %%
# prefix so these IDs are unique from other datasets
infection_df["ID"] = "VariantInfection_" + infection_df["ID"].astype(str)

# %%
infection_df.groupby(["vaccinated", "vaccine"]).size()

# %%
infection_df["mAb_cocktail"].fillna("no", inplace=True)

# %%
# Remove mAb-treated patients - the external antibodies may interfere with our measurements
infection_df = infection_df[
    (infection_df["vaccinated"] != "mAb") & (infection_df["mAb_cocktail"] != "yes")
].copy()

# %%
assert not infection_df["vaccinated"].isna().any()
assert not infection_df["vaccine"].isna().any()
assert not infection_df["mAb_cocktail"].isna().any()

# %%
infection_df.groupby(["vaccinated", "vaccine", "mAb_cocktail"]).size()

# %%

# %%
infection_df.head()

# %%
infection_df.columns

# %%
infection_df["variant"].value_counts()

# %%
infection_df["variant"] = infection_df["variant"].str.title()

# %%
infection_df["variant"].value_counts()

# %%

# %%
infection_df["vaccine"].value_counts()

# %%
infection_df["vaccine"] = infection_df["vaccine"].replace(
    {"Pfizer": "mRNA", "Moderna": "mRNA"}
)

# %%
infection_df["vaccine"].value_counts()

# %%

# %%
measurement_cols = infection_df.columns[infection_df.columns.str.contains("AU")]
helpers.confirm_all_measurement_columns_are_present(measurement_cols)
measurement_cols

# %%

# %%
# only one sample per patient
assert all(infection_df["ID"].value_counts() == 1)

# %%

# %% [markdown]
# # reformat

# %%

# %%
# Separate into groups by variant and by vaccination status if any

infection_df_obs = infection_df[
    [
        "ID",
        "variant",
        "vaccine",
    ]
].copy()

infection_df_obs["Status"] = "Variant Infection" + " - " + infection_df_obs["variant"]
# add vaccine info suffix if vaccinated
infection_df_obs.loc[infection_df_obs["vaccine"] != "Not vaccinated", "Status"] += (
    " - "
    + infection_df_obs.loc[infection_df_obs["vaccine"] != "Not vaccinated", "vaccine"]
    + " vaccinated"
)

# capitalize columns for consistency - only the first letter
# and set index
infection_df_obs = infection_df_obs.rename(
    columns=lambda s: s[0].upper() + s[1:]
).set_index("ID")

# anndata wants string index
infection_df_obs.index = infection_df_obs.index.astype(str)

infection_df_obs

# %%

# %%
# extract AU measurement cols, and set index to match obs
infection_df_X = (
    infection_df[measurement_cols]
    .rename(columns=lambda col: col.replace("_AU", ""))
    .set_index(infection_df_obs.index)
)

infection_df_X

# %%

# %%

# %%
adata = anndata.AnnData(X=infection_df_X, obs=infection_df_obs)
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
adata.var["timepoint"] = "day 28 / week 4"
adata.var

# %%
# create combined name
adata.var["combined_name"] = adata.var.apply("_".join, axis=1)
adata.var = adata.var.set_index("combined_name")

# %%
adata.var

# %%

# %%

# %%
adata

# %%
# filter down variant infection types
adata.obs["Status"].value_counts()

# %%
# filter down variant infection types
adata = adata[
    adata.obs["Status"].isin(
        [
            "Variant Infection - Delta",
            "Variant Infection - Delta - mRNA vaccinated",
            "Variant Infection - Alpha",
        ]
    )
]

# %%
adata

# %%

# %%

# %%
adata.write(f"{config.paths.generated_data_dir}/partial.variant_infections.h5")

# %%

# %%

# %%
