import seaborn as sns
from genetools.palette import HueValueStyle

# Confirm plate 11 measurement column names
expected_measurement_columns = [
    "Epsilon_AU",
    "Beta_AU",
    "Iota_AU",
    "Gamma_AU",
    "B.1.526.2_AU",
    "Alpha_AU",
    "P.3_AU",
    "Kappa_AU",
    "Delta_AU",
    "Wuhan_AU",
]


def confirm_all_measurement_columns_are_present(measurement_cols):
    diff = set(measurement_cols).symmetric_difference(set(expected_measurement_columns))
    if len(diff) != 0:
        raise ValueError(f"Wrong measurement column set: difference = {diff}")
    return True


# Define plotting palette
# sRGB values grabbed from Katharina's figures
palette_dict = {
    "Outpt / Admit": (249, 141, 60),
    "Infection": (249, 141, 60),
    "Wuhan Infection - Admit": (249, 141, 60),
    "Wuhan Infection - Moderate": (249, 141, 60),
    "ICU / Death": (203, 29, 29),
}

# convert to float
for key, color_tuple in palette_dict.items():
    palette_dict[key] = tuple(val / 255 for val in color_tuple)

# add other colors
# map certain hue values to different marker styles, with size scale factors
# for unfilled shapes, set facecolors="none" but set edgecolors to the main color.
default_palette = sns.color_palette("bright")
palette_dict.update(
    {
        "Control": default_palette[7],
        "Wuhan Infection - Outpatient": default_palette[2],
        "Infection": HueValueStyle(
            color=palette_dict["Infection"], zorder=20, marker_size_scale_factor=1.5
        ),
        "Wuhan Infection - Mild": default_palette[2],
        "Wuhan Infection - ICU": HueValueStyle(color=default_palette[3], zorder=3),
        "Wuhan Infection - Severe": HueValueStyle(color=default_palette[3], zorder=3),
        "Pfizer-Pfizer (Stanford)": HueValueStyle(
            color=default_palette[0],
            edgecolors=default_palette[0],
            facecolors="none",
            marker="^",
            marker_size_scale_factor=1.5,
            linewidths=1.5,
        ),
        "Pfizer-Pfizer (Stanford), CoV2+": HueValueStyle(
            color=default_palette[9],
            zorder=5,
            marker="*",
            marker_size_scale_factor=4.5,
            legend_size_scale_factor=2.0,
            linewidths=1.5,
            alpha=1.0,
        ),
        "Sinopharm-Sinopharm": HueValueStyle(
            color=default_palette[5],
            edgecolors=default_palette[5],
            facecolors="none",
            marker="s",
            marker_size_scale_factor=1.5,
            linewidths=2.0,
            zorder=10,
        ),
        "Sputnik V-Sputnik V": HueValueStyle(
            color=default_palette[6],
            edgecolors=default_palette[6],
            facecolors="none",
            marker="s",
            marker_size_scale_factor=1.5,
            linewidths=2.0,
            zorder=10,
        ),
        "AstraZeneca-AstraZeneca": HueValueStyle(
            color=sns.color_palette("viridis")[2],
            edgecolors=sns.color_palette("viridis")[2],
            facecolors="none",
            marker="s",
            marker_size_scale_factor=1.5,
            linewidths=2.0,
            zorder=10,
        ),
        "Pfizer-Pfizer (Mongolia)": HueValueStyle(
            color=default_palette[9],
            edgecolors=default_palette[9],
            facecolors="none",
            marker="^",
            marker_size_scale_factor=1.5,
            linewidths=1.5,
            zorder=10,
        ),
        "Variant Infection - Alpha": HueValueStyle(
            color=default_palette[6],
            edgecolors=default_palette[6],
            facecolors="none",
            marker="s",
            marker_size_scale_factor=1.5,
            linewidths=2.0,
            zorder=12,
        ),
        "Variant Infection - Delta": HueValueStyle(
            color=sns.color_palette("viridis")[2],
            edgecolors=sns.color_palette("viridis")[2],
            facecolors="none",
            marker="s",
            marker_size_scale_factor=1.5,
            linewidths=2.0,
            zorder=12,
        ),
    }
)
