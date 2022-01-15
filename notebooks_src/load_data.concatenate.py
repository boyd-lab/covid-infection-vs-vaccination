# -*- coding: utf-8 -*-
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

# %% [markdown]
# # Load data

# %%
adata_sources = {
    "Pfizer_vaccine_coronavirus_plate": f"{config.paths.generated_data_dir}/partial.pfizer_vaccine.coronavirus_plate.h5",
    "Pfizer_vaccine_variant_plate": f"{config.paths.generated_data_dir}/partial.pfizer_vaccine.variant_plate.h5",
    "infection_cohort2_coronavirus_plate": f"{config.paths.generated_data_dir}/partial.infection_cohort2.coronavirus_plate.h5",
    "infection_cohort2_variant_plate": f"{config.paths.generated_data_dir}/partial.infection_cohort2.variant_plate.h5",
    "Mongolia_vaccine_variant_plate": f"{config.paths.generated_data_dir}/partial.mongolia_vaccines.h5",
    "infection_cohort1": f"{config.paths.generated_data_dir}/partial.infection_cohort1.h5",
    "variant_infections": f"{config.paths.generated_data_dir}/partial.variant_infections.h5",
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

# %% [markdown]
# ## Combine coronavirus and variant plate anndatas from same cohort

# %%
adatas.keys()

# %%
adatas_by_cohort = {
    "Mongolia_vaccine_variant_plate": adatas["Mongolia_vaccine_variant_plate"],
    "infection_cohort1": adatas["infection_cohort1"],
    "variant_infections": adatas["variant_infections"],
}

# %%
adatas_by_cohort["Pfizer_all"] = anndata.concat(
    [
        adatas["Pfizer_vaccine_coronavirus_plate"],
        adatas["Pfizer_vaccine_variant_plate"],
    ],
    join="outer",
    merge="first",
    axis=1,
)
adatas["Pfizer_vaccine_coronavirus_plate"], adatas[
    "Pfizer_vaccine_variant_plate"
], adatas_by_cohort["Pfizer_all"]

# %%
adatas_by_cohort["infection_cohort2_all"] = anndata.concat(
    [
        adatas["infection_cohort2_coronavirus_plate"],
        adatas["infection_cohort2_variant_plate"],
    ],
    join="outer",
    merge="first",
    axis=1,
)
adatas["infection_cohort2_coronavirus_plate"], adatas[
    "infection_cohort2_variant_plate"
], adatas_by_cohort["infection_cohort2_all"]

# %%
adatas_by_cohort

# %%
for name, adata in adatas_by_cohort.items():
    display(Markdown(f"## {name}"))
    display(adata.var)

# %%

# %%
# Confirm no obs names overlap between datasets
import itertools

for ((name_a, adata_a), (name_b, adata_b)) in itertools.combinations(
    adatas_by_cohort.items(), 2
):
    intersection_of_obsnames = set.intersection(
        set(adata_a.obs_names), set(adata_b.obs_names)
    )
    print(name_a, name_b, intersection_of_obsnames)
    assert len(intersection_of_obsnames) == 0


# %%
# Note issue: https://github.com/theislab/anndata/issues/614
# This doesn't work - adata.var has a lot of NaNs.

# adata_full = anndata.concat(
#     adatas_by_cohort.values(),
#     join="outer",
#     merge="first",
#     axis=0,
#     label="source_cohort",
#     keys=adatas_by_cohort.keys(),
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
for name, adata in adatas_by_cohort.items():
    adata.obs["source_cohort"] = name

# now merge
adata_full = merge_anndatas(
    adatas_by_cohort.values(),
    var_col_join=["virus", "target", "variant_plate_type", "antibody", "timepoint"],
)

adata_full

# %%
adata_full.var

# %%
adata_full.obs

# %%

# %%

# %%
adata_full.obs.columns

# %%
assert not adata_full.obs_names.duplicated().any()

# %%
adata_full.obs[adata_full.obs["Timepoint"].isna()]["source_cohort"].value_counts()


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
    {"Asymptomatic": "Mild", "Critical": "Severe"}
)
adata_full.obs["Exposure"].value_counts()

# %%
adata_full.obs["Exposure"] = adata_full.obs["Exposure"].replace(
    {
        "Mild": "Wuhan Infection - Mild",
        "Moderate": "Wuhan Infection - Moderate",
        "Severe": "Wuhan Infection - Severe",
    }
)
adata_full.obs["Exposure"].value_counts()

# %%
adata_full.obs["Exposure"] = adata_full.obs["Exposure"].replace(
    {
        "Pfizer-Pfizer": "Pfizer-Pfizer (Mongolia)",
        "Vaccinee": "Pfizer-Pfizer (Stanford)",
        "Vaccinee (CoV2+)": "Pfizer-Pfizer (Stanford), CoV2+",
        "AZ-AZ": "AstraZeneca-AstraZeneca",
    }
)
adata_full.obs["Exposure"].value_counts()

# %%

# %%

# %%
patient_types = [
    "Admit",
    "Wuhan Infection - Mild",
    "ICU",
    "Outpatient",
    "Variant Infection - Delta",
    "Variant Infection - Alpha",
    "Wuhan Infection - Moderate",
    "Wuhan Infection - Severe",
]
adata_full.obs["Exposure Type"] = adata_full.obs["Exposure"].replace(
    {k: "Infection" for k in patient_types}
)
adata_full.obs["Exposure Type"].value_counts()

# %%

# %%


# %%

# %% [markdown]
# # Variant vs Wuhan ratios
#
# - Variant ratios calculated as Wuhan divided by variant value, from plate 11 data.
# - Variant measurements checked against a variant plate cutoff threshold determined by the sample's exposure type; e.g. Delta infeections had their Delta measurement checked against the Delta cutoff.

# %%
adata_full.var["combined_name"] = adata_full.var_names
adata_full.var

# %%
igg_rbd_level_columns = adata_full.var.query(
    'target == "RBD" & antibody == "IgG" & variant_plate_type == "Variant"'
)
igg_rbd_level_columns

# %%
igg_rbd_level_columns["virus"].value_counts()

# %%
orig_data = adata_full.to_df()
ratios = []
ratios_annot = []

for timepoint, igg_rbd_level_columns_by_timepoint in igg_rbd_level_columns.groupby(
    "timepoint"
):
    # Compute ratios for each timepoint

    # Numerator is always the wuhan level
    numerator_column = (
        igg_rbd_level_columns_by_timepoint[
            igg_rbd_level_columns_by_timepoint["virus"] == "Wuhan"
        ]
        .squeeze()
        .name
    )

    # Many denominators: each non-wuhan measurement
    denominator_columns = (
        igg_rbd_level_columns_by_timepoint[
            igg_rbd_level_columns_by_timepoint.index != numerator_column
        ]
        .reset_index()
        .rename(columns={"varname": "denominator"})
        .assign(numerator=numerator_column)
    )

    # Make combined name
    denominator_columns["combined_name"] = (
        "ratio:"
        + denominator_columns["numerator"]
        + ":"
        + denominator_columns["denominator"]
    )

    print(timepoint, numerator_column, denominator_columns["denominator"].values)

    # Create one numerator/denominator entry for each denominator
    for _, denominator_columns_annot in denominator_columns.iterrows():
        denominator_column = denominator_columns_annot["denominator"]
        ratios.append(orig_data[numerator_column] / orig_data[denominator_column])
        ratios_annot.append(denominator_columns_annot)

# %%
# Combine all ratios across all timepoints (blanks if sample not measured in a particular timepoint)
ratios = pd.DataFrame(ratios).T
ratios

# %%
# Obs match
assert ratios.index.equals(adata_full.obs_names)

# %%
# Now label the columns with the (numerator, denominator, timepoint) info
ratios_annot = pd.DataFrame(ratios_annot).set_index("combined_name")
ratios_annot

# %%

# %%

# %%
# Stitch those annotations into a ratios anndata
ratios_adata = anndata.AnnData(ratios, var=ratios_annot, obs=adata_full.obs.copy())
ratios_adata

# %%
ratios_adata.var

# %%
# here's a preview of what we created
ratios_adata.to_df()

# %%

# %% [markdown]
# # Add cutoffs for variant ratios

# %% [markdown]
# ## Apply cutoffs
#
# We have:
#
# - Patients exposed to variants or original Wuhan virus. We know whether the patient was exposed to Delta, or Alpha, or Wuhan.
# - Vaccinees -- all exposed to Wuhan type virus.
#
# Compare only the relevant column for each patient. For example, for a Delta-exposed patient, look for whether they are violating the Delta cutoff.

# %%
plate11_cutoffs = (
    pd.read_csv(f"{config.paths.data_dir}/cutoff_p11.csv")
    .rename(columns=lambda s: s.replace("_AU", ""))
    .T.rename(columns={0: "P11_cutoff"})
)
plate11_cutoffs

# %%
adata_full.var[adata_full.var["variant_plate_type"] == "Variant"]

# %%

# %%
# merge on measurement_original_column_name to find column name in adata
plate11_cutoffs_df = pd.merge(
    adata_full.var[adata_full.var["variant_plate_type"] == "Variant"],
    plate11_cutoffs,
    left_on="virus",
    right_index=True,
    how="inner",
    validate="m:1",
)
plate11_cutoffs_df

# %%
plate11_cutoffs_df["virus"].unique()

# %%
# figure out which column we want to check for each patient
adata_full.obs["Exposure"].value_counts()

# %%
# figure out which measurement cutoff we want to check for each patient

# We have:
# Patients exposed to variants or Wuhan -> know whether they’re exposed to Delta, or Alpha, or Wuhan
# Vaccinees -> exposed to Wuhan

# Compare only the relevant column for each patient. E.g. if Delta-exposed patient, look for whether they are violating the Delta cutoff.

adata_full.obs["cutoff_column_to_check"] = adata_full.obs["Exposure"].replace(
    {
        "Wuhan Infection - Mild": "Wuhan",
        "Pfizer-Pfizer (Stanford)": "Wuhan",
        "AstraZeneca-AstraZeneca": "Wuhan",
        "Admit": "Wuhan",
        "Sinopharm-Sinopharm": "Wuhan",
        "ICU": "Wuhan",
        "Sputnik V-Sputnik V": "Wuhan",
        "Pfizer-Pfizer (Mongolia)": "Wuhan",
        "Outpatient": "Wuhan",
        "Wuhan Infection - Severe": "Wuhan",
        "Wuhan Infection - Moderate": "Wuhan",
        "Pfizer-Pfizer (Stanford), CoV2+": "Wuhan",
        "Control": "Wuhan",
        "Variant Infection - Delta": "Delta",
        "Variant Infection - Delta - mRNA vaccinated": "Delta",
        "Variant Infection - Alpha": "Alpha",
    }
)
assert all(
    col in plate11_cutoffs_df["virus"].unique()
    for col in adata_full.obs["cutoff_column_to_check"]
)
adata_full.obs["cutoff_column_to_check"].value_counts()

# %%
# For each row (each sample), mark which cutoffs are violated - we do this one cutoff at a time:

is_below_cutoff = {}
for varname, cutoff_info in plate11_cutoffs_df.iterrows():
    # For each cutoff, find rows where it is violated
    cutoff_value = cutoff_info["P11_cutoff"]
    varname_virus = cutoff_info["virus"]
    print(varname, varname_virus, cutoff_value)

    # Check only rows (samples) that match the current column. E.g. if varname is delta (at some timepoint), filter to Delta patients only
    is_below_cutoff[varname] = (
        adata_full[adata_full.obs["cutoff_column_to_check"] == varname_virus].to_df()[
            varname
        ]
        < cutoff_value
    )

    # indicate that all other rows do not violate this cutoff
    is_below_cutoff[varname] = (
        is_below_cutoff[varname].reindex(adata_full.obs_names).fillna(False)
    )

# combine data from all cutoffs
is_below_cutoff = pd.DataFrame(is_below_cutoff)

# confirm indexed the same as adata_full
assert np.array_equal(is_below_cutoff.index, adata_full.obs_names)
is_below_cutoff

# %%
# how many cutoffs are violated per patient?
is_below_cutoff.astype(int).sum(axis=1).value_counts()

# %%

# %%
# get whether each patient has violated any cutoff
any_cutoffs_violated = is_below_cutoff.any(axis=1)

# most patients have not violated a cutoff at any time point:
any_cutoffs_violated.value_counts()

# %%
# set obs
adata_full.obs["any_p11_cutoffs_violated"] = any_cutoffs_violated
# update other obs too, which is just a copy
ratios_adata.obs = adata_full.obs.copy()

adata_full.obs.head()

# %%
# who is violating cutoffs?
adata_full.obs[adata_full.obs["any_p11_cutoffs_violated"]]["Status"].value_counts()

# %%

# %%
# within each timepoint: which patients have violated any field's cutoff?
for timepoint, plate11_cutoffs_df_for_timepoint in plate11_cutoffs_df.groupby(
    "timepoint"
):
    # who has or hasn't have violated a cutoff at this time point:
    any_cutoffs_violated_at_timepoint = is_below_cutoff[
        plate11_cutoffs_df_for_timepoint.index
    ].any(axis=1)

    print(timepoint, any_cutoffs_violated_at_timepoint.value_counts())

    # set obs
    adata_full.obs[
        f"any_p11_cutoffs_violated:{timepoint}"
    ] = any_cutoffs_violated_at_timepoint

# %%
# update other obs too, which is just a copy
ratios_adata.obs = adata_full.obs.copy()

# %%
# who is violating cutoffs at each timepoint?
for timepoint, _ in plate11_cutoffs_df.groupby("timepoint"):
    print(timepoint)
    print(
        adata_full.obs[adata_full.obs[f"any_p11_cutoffs_violated:{timepoint}"]][
            "Status"
        ].value_counts()
    )
    print("---")
    print()

# %%

# %%

# %%

# %% [markdown]
# # Concatenate

# %%
adata_full, ratios_adata

# %%
assert adata_full.obs_names.equals(ratios_adata.obs_names)

# %%
adata_merge = anndata.concat(
    [adata_full, ratios_adata], axis=1, join="outer", merge="first"
)
adata_merge

# %%
adata_merge.var

# %%
adata_merge.obs

# %%
assert adata_merge.obs.equals(adata_full.obs)
assert adata_merge.obs.equals(ratios_adata.obs)

# %%
# confirm no nans
assert not adata_merge.obs["any_p11_cutoffs_violated"].isna().any()

# %%
cutoff_violation_cols = adata_merge.obs.columns[
    adata_merge.obs.columns.str.startswith("any_p11_cutoffs_violated")
]
print(cutoff_violation_cols)
for col in cutoff_violation_cols:
    assert not adata_merge.obs[col].isna().any(), col

# %%

# %%

# %% [markdown]
# # Export

# %%
adata_full.write(f"{config.paths.generated_data_dir}/combined_data.h5")

# %%
ratios_adata.write(f"{config.paths.generated_data_dir}/variant_ratios.h5")

# %%
adata_merge.write(
    f"{config.paths.generated_data_dir}/wuhan_measurements_and_variant_ratios.h5"
)

# %%

# %%

# %%
