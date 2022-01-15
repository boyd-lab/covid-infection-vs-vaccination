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
from numpy.testing import assert_array_equal

# %%

# %% [markdown]
# # Load data

# %%
infection_df = pd.read_csv(f"{config.paths.data_dir}/infection_cohort1.csv")
infection_df

# %%

# %%

# %% [markdown]
# # infection data reshape

# %%
infection_df

# %%
infection_df[infection_df["Status"].isna()]

# %%
infection_df.dropna(how="any", subset=["Patient_ID", "Status"], inplace=True)
infection_df

# %%

# %%
melt_id_vars = [
    "sample_id",
    "plate_id",
    "Patient_ID",
    "Status",
    "Death",
    "Days_since_symptom_onset",
]

# %%
infection_df = pd.melt(
    infection_df,
    id_vars=melt_id_vars,
    var_name="measurement",
)
infection_df

# %%
infection_df.dtypes

# %%
infection_df["measurement"].unique()

# %%
# there are some non-numeric values in "value"
infection_df["value"].apply(np.isreal).value_counts()

# %%
"#VALUE!" in infection_df["value"].values, "#DIV/0!" in infection_df["value"].values

# %%
# convert to float and switch non-numeric values to nan
infection_df["value"] = pd.to_numeric(infection_df["value"], errors="coerce")
infection_df["value"].isna().value_counts()

# %%

# %%
infection_df["plate_id"].value_counts()

# %%

# %%
infection_df = infection_df[
    (infection_df["measurement"].str.contains("IgG"))
    | (infection_df["measurement"].str.contains("IgA"))
    | (infection_df["measurement"].str.contains("IgM"))
]
infection_df

# %%
infection_df["measurement"].unique()

# %%
infection_df = infection_df[
    infection_df["measurement"].str.lower().str.contains("mean")
].copy()
infection_df

# %%
infection_df["measurement"].value_counts()

# %%
infection_df["measurement_original_column_name"] = infection_df["measurement"].copy()

# %%
infection_df["measurement"] = (
    infection_df["measurement"].str.replace("_Mean", "").str.replace("_mean", "")
)
infection_df["measurement"].value_counts()

# %%
# extract measurement parts
measurement_parts = (
    infection_df["measurement"]
    .str.split("_", expand=True)
    .rename(columns={0: "virus", 1: "target", 2: "antibody"})
    .assign(variant_plate_type="Wuhan")
    .apply(lambda col: col.str.strip())
)
measurement_parts

# %%
measurement_parts["virus"] = measurement_parts["virus"].replace({"CoV2": "Wuhan"})

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
# filter out ELISA readings
infection_df = infection_df[infection_df["virus"] != "ELISA"]
infection_df["virus"].value_counts()

# %%
infection_df["plate_id"].value_counts()

# %%
infection_df["Days_since_symptom_onset"].value_counts()

# %%
sorted(infection_df["Days_since_symptom_onset"].unique())

# %%

# %%

# %%
# extract patient status - take first row for each patient
infection_patient_status = infection_df.groupby("Patient_ID").first()[
    ["Status", "Death"]
]
infection_patient_status

# %%

# %%
infection_patient_status.drop_duplicates()

# %%
# change "admit? Unclear if was in ICU" to "admit". And capitalize them all
infection_patient_status["Status"] = infection_patient_status["Status"].replace(
    {"admit": "Admit", "icu": "ICU", "admit? Unclear if was in ICU": "Admit"}
)

# %%
infection_patient_status["Status"].value_counts()

# %%
# doesn't include death NaNs
infection_patient_status.groupby(["Status", "Death"]).size()

# %%
# what are the death NaNs: variant infections and prepandemic control
infection_patient_status[infection_patient_status["Death"].isna()][
    "Status"
].value_counts()

# %%
infection_patient_status[infection_patient_status["Status"] == "Control"]

# %%

# %%
print(infection_patient_status.shape)
infection_patient_status.head()

# %%
infection_patient_status["Status"].value_counts()

# %%

# %%
# patients are distributed across plates - almost all patients only are analyzed on a single plate
infection_df.groupby("Patient_ID")["plate_id"].nunique().value_counts()


# %%

# %%
def process_infection_timepoint(df_partial, timepoint):
    # as weeks
    (df_partial["Days_since_symptom_onset"] / 7).hist()
    plt.xlabel("weeks")

    # measurements of patient-virus-target-antibody combos are repeated across days - some patients have multiple samples in this time period
    # but not repeated within a day
    # i.e. for each day, one patient-virus-target-antibody measurement made — and only on one single plate
    assert all(
        df_partial.groupby(
            [
                "Patient_ID",
                "virus",
                "target",
                "antibody",
                "Days_since_symptom_onset",
                "variant_plate_type",
            ]
        ).size()
        == 1
    )

    # for each individual, take mean of samples in this time frame (across plates)
    # (this combines measurements across plates)
    df_partial = (
        df_partial.groupby(
            [
                "Patient_ID",
                "virus",
                "target",
                "variant_plate_type",
                "antibody",
                "measurement_original_column_name",
            ],
            sort=False,
        )["value"]
        .mean()
        .reset_index()
    )

    # sanity check
    assert all(
        df_partial.groupby(
            ["Patient_ID", "virus", "target", "antibody", "variant_plate_type"]
        ).size()
        == 1
    )

    # pivot to matrix
    infection_df_pivot = pd.pivot(
        df_partial,
        index="Patient_ID",
        columns=[
            "virus",
            "target",
            "variant_plate_type",
            "antibody",
            "measurement_original_column_name",
        ],
        values="value",
    )

    #     # drop patients with any NaNs
    #     infection_df_pivot = infection_df_pivot.dropna(how="any")
    #     assert not infection_df_pivot.isna().any().any()

    ## set column names
    variable_info = infection_df_pivot.columns.to_frame().reset_index(drop=True)
    # create combined name
    variable_info["timepoint"] = timepoint
    variable_info["combined_name"] = variable_info[
        variable_info.columns.drop("measurement_original_column_name")
    ].apply("_".join, axis=1)
    variable_info = variable_info.set_index("combined_name")

    # set var names
    infection_df_pivot.columns = variable_info.index.copy()

    return infection_df_pivot, variable_info


# %%

# %%
infection_df["Days_since_symptom_onset"].hist()

# %%
assert (
    infection_df[infection_df["Status"] == "Control"]["Days_since_symptom_onset"]
    .isna()
    .all()
)
assert (
    infection_df[infection_df["Status"].str.startswith("Variant infection")][
        "Days_since_symptom_onset"
    ]
    .isna()
    .all()
)

# %%
# map infection timepoints to global timepoint label
map_infection_to_global_timepoint_labels = {
    # infection data points with Days_since_symptom_onset N/A can be control samples
    "day 0 / pre-pandemic": infection_df[
        (infection_df["Status"] == "Control")
        & (infection_df["Days_since_symptom_onset"].isna())
    ],
    # exclude infection data points with Days_since_symptom_onset<= 0
    "day 7 / week 1": infection_df[
        (infection_df["Days_since_symptom_onset"] >= 0)
        & (infection_df["Days_since_symptom_onset"] <= 1 * 7)
    ],
    "day 21 / weeks 2&3": infection_df[
        (infection_df["Days_since_symptom_onset"] > 1 * 7)
        & (infection_df["Days_since_symptom_onset"] <= 3 * 7)
    ],
    "day 28 / week 4": infection_df[
        (infection_df["Days_since_symptom_onset"] > 3 * 7)
        & (infection_df["Days_since_symptom_onset"] <= 4 * 7)
    ],
    "day 42 / weeks 5&6": infection_df[
        (infection_df["Days_since_symptom_onset"] > 4 * 7)
        & (infection_df["Days_since_symptom_onset"] <= 6 * 7)
    ],
    "week 7 and later / 3 months": infection_df[
        infection_df["Days_since_symptom_onset"] > 6 * 7
    ],
}

for (
    associated_global_timepoint_label,
    df_partial,
) in map_infection_to_global_timepoint_labels.items():
    print(df_partial.shape, associated_global_timepoint_label)

# %%
# sanity check that every infection datapoint with a non-nan and >=0 timepoint is accounted for in these splits
assert (
    sum(
        [
            df_partial.shape[0]
            for df_partial in map_infection_to_global_timepoint_labels.values()
        ]
    )
    + infection_df[infection_df["Status"] != "Control"]["Days_since_symptom_onset"]
    .isna()
    .sum()
    + (
        infection_df[infection_df["Status"] != "Control"]["Days_since_symptom_onset"]
        < 0
    ).sum()
    == infection_df.shape[0]
)

# %%
X_partial = []
var_partial = []
for (
    associated_global_timepoint_label,
    df_partial,
) in map_infection_to_global_timepoint_labels.items():
    print(associated_global_timepoint_label)
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
# reorder participants info to match
infection_patient_status = infection_patient_status.loc[infection_df_pivot.index]
infection_patient_status

# %%
infection_patient_status.columns

# %%

# %%
adata_infection = anndata.AnnData(
    X=infection_df_pivot, obs=infection_patient_status, var=variable_info
)
adata_infection

# %%
adata_infection.var

# %%
adata_infection.obs["Status"].value_counts()

# %%

# %%
# Remove "Control"

# %%
adata_infection = adata_infection[adata_infection.obs["Status"] != "Control"]

# %%
adata_infection.obs["Status"].value_counts()

# %%

# %%
adata_infection.write(f"{config.paths.generated_data_dir}/partial.infection_cohort1.h5")

# %%
