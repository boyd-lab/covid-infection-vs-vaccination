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

# %%
from IPython.display import display, Markdown
from numpy.testing import assert_array_equal

# %%

# %%

# %%

# %% [markdown]
# # Concatenate a subset of datasets: coronavirus plate only, for `infection_cohort1` patients and for Pfizer vaccinees (original data that includes IgM and IgA -- new only includes IgG)
#
# # Load data

# %%
adata_sources = {
    "Pfizer_vaccine": f"{config.paths.generated_data_dir}/partial.pfizer_vaccine.coronavirus_plate.original.h5",
    "infection_cohort1": f"{config.paths.generated_data_dir}/partial.infection_cohort1.h5",
}

# %%
adatas = {key: sc.read(val) for key, val in adata_sources.items()}
adatas

# %%

# %% [markdown]
# # combine datasets and align columns

# %%
for name, adata in adatas.items():
    display(Markdown(f"## {name}"))
    display(adata.var)
    print(adata.var["timepoint"].unique().tolist())

# %%

# %%
# Confirm no obs names overlap between datasets
import itertools

for ((name_a, adata_a), (name_b, adata_b)) in itertools.combinations(adatas.items(), 2):
    intersection_of_obsnames = set.intersection(
        set(adata_a.obs_names), set(adata_b.obs_names)
    )
    print(name_a, name_b, intersection_of_obsnames)
    assert len(intersection_of_obsnames) == 0


# %%
# # Note issue: https://github.com/theislab/anndata/issues/614
# # This doesn't work - adata.var has a lot of NaNs.

# adata_full = anndata.concat(
#     adatas.values(),
#     join="outer",
#     merge="first",
#     axis=0,
#     label="source_cohort",
#     keys=adatas.keys(),
# )
# adata_full

# %%
def _merge_two_anndatas(ad1, ad2, var_col_join):
    """concatenate two anndatas. they must have different obsnames.
    some of their vars can intersect, and will be combined.
    we will concat along var, specifically using only the [var_col_join] columns to describe each variable, as well as the var_name
    """
    # confirm obsnames are distinct
    if len(set(ad1.obs_names).intersection(ad2.obs_names)) > 0:
        raise ValueError("Obsnames intersect")

    # we will concat along var, specifically using only the [var_col_join] columns to describe each variable, as well as the var_name
    if "varname" in var_col_join:
        # TODO: relax this
        raise ValueError("Cannot use varname as a var col - will be overwritten")

    def _get_df_from_anndata(adata, var_cols):
        df = adata.to_df()
        df.columns = adata.var[var_cols].assign(varname=adata.var_names)
        return df

    df1 = _get_df_from_anndata(ad1, var_col_join)
    df2 = _get_df_from_anndata(ad2, var_col_join)

    df_concat = pd.concat([df1, df2], axis=0)
    if df_concat.shape[0] != df1.shape[0] + df2.shape[0]:
        raise ValueError("Concat produced unexpected number of rows")

    new_var = df_concat.columns.to_frame(index=False)
    new_var.columns = var_col_join + ["varname"]
    # recover varname
    new_var = new_var.set_index("varname")
    df_concat.columns = new_var.index

    new_obs = pd.concat([ad1.obs, ad2.obs], axis=0)
    if not np.array_equal(new_obs.index, df_concat.index):
        raise ValueError("Concat unexpectedly rearranged rows")

    if new_var.duplicated().any():
        raise ValueError("Some var rows are duplicated - unexpected")

    return anndata.AnnData(df_concat, var=new_var, obs=new_obs)


from functools import reduce


def merge_anndatas(
    adatas,
    var_col_join,
):
    """progresively merge a list of anndatas"""
    return reduce(
        lambda x, y: _merge_two_anndatas(x, y, var_col_join=var_col_join), adatas
    )


# %%
# use our workaround
# first, label each adata with a source cohort key
for name, adata in adatas.items():
    adata.obs["source_cohort"] = name

# now merge
adata_full = merge_anndatas(
    adatas.values(),
    var_col_join=["virus", "target", "variant_plate_type", "antibody", "timepoint"],
)

adata_full

# %%
adata_full.var

# %%
adata_full.obs

# %%

# %%
adatas.keys()

# %%
var_names_in_common = adatas["infection_cohort1"].var_names.intersection(
    adatas["Pfizer_vaccine"].var_names
)
len(var_names_in_common), var_names_in_common

# %%
var_names_different = adatas["infection_cohort1"].var_names.difference(
    adatas["Pfizer_vaccine"].var_names
)
len(var_names_different), var_names_different

# %%
adata_full

# %%
# subset to columns in common
adata_full = adata_full[:, var_names_in_common].copy()
adata_full

# %%

# %%

# %%

# %%
assert not adata_full.obs_names.duplicated().any()

# %%
assert not adata_full.obs["Status"].isna().any()

# %%

# %%
adata_full.var["timepoint"].unique().tolist()

# %%

# %%

# %% [markdown]
# # Expand granularity of Status obs column, and add any other hue columns

# %%
adata_full.obs["Status"].value_counts()

# %%
adata_full.obs["Exposure"] = adata_full.obs["Status"].copy()

# %%
adata_full.obs.loc[
    (adata_full.obs["Exposure"] == "Vaccinee")
    & (adata_full.obs["COVID Positive Ever?"] != "No")
    & ~(adata_full.obs["COVID Positive Ever?"].isna()),
    "Exposure",
] = "Vaccinee (CoV2+)"
adata_full.obs["Exposure"].value_counts()

# %%
adata_full.obs["Exposure"] = adata_full.obs["Exposure"].replace(
    {
        "Admit": "Wuhan Infection - Admit",
        "ICU": "Wuhan Infection - ICU",
        "Outpatient": "Wuhan Infection - Outpatient",
    }
)
adata_full.obs["Exposure"].value_counts()

# %%
adata_full.obs["Exposure"] = adata_full.obs["Exposure"].replace(
    {
        "Vaccinee": "Pfizer-Pfizer (Stanford)",
        "Vaccinee (CoV2+)": "Pfizer-Pfizer (Stanford), CoV2+",
    }
)
adata_full.obs["Exposure"].value_counts()

# %%

# %%

# %%
patient_types = [
    "Wuhan Infection - Admit",
    "Wuhan Infection - ICU",
    "Wuhan Infection - Outpatient",
]
adata_full.obs["Exposure Type"] = adata_full.obs["Exposure"].replace(
    {k: "Infection" for k in patient_types}
)
adata_full.obs["Exposure Type"].value_counts()


# %%

# %% [markdown]
# # Export

# %%
adata_full.write(f"{config.paths.generated_data_dir}/coronavirus_plate_only.subset.h5")

# %%
