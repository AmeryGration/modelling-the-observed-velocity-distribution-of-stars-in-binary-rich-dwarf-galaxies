#!/usr/bin/env python

"""Compute binary kinematics"""

import warnings
import yaml
import numpy as np
import pandas as pd
import dyad

########################################################################
# Load variables from configuration file
########################################################################
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

population = config["population"]
data_dir = population + "/data/"

########################################################################
# Load binary population data
########################################################################
df_properties = pd.read_csv(data_dir + "binary_properties.csv")
# df_properties = df_properties.sample(1_000_000, ignore_index=True)
try:
    weight = df_properties["weight"]
except KeyError:
    warnings.warn(
        "data frame does not contain a weight column.", UserWarning
    )
    df_properties["weight"] = 1.
    weight = df_properties["weight"]

########################################################################
# Compute binary kinematics and save to file
########################################################################
binary = dyad.TwoBody(
    df_properties["mass_1/Msun"],
    df_properties["mass_ratio"],
    (
        df_properties["semimajor_axis/AU"]
        *df_properties["mass_ratio"]
        /(1. + df_properties["mass_ratio"])
    ),
    df_properties["eccentricity"],
    df_properties["longitude_of_ascending_node"],
    df_properties["inclination"],
    df_properties["argument_of_pericentre"],
)

cols = [
    "mass_1/Msun",
    "x_1/AU",
    "y_1/AU",
    "z_1/AU",
    "v_1x/km/s",
    "v_1y/km/s",
    "v_1z/km/s",
    "radius_1/AU",
    "v_1/km/s",
    "mass_2/Msun",
    "x_2/AU",
    "y_2/AU",
    "z_2/AU",
    "v_2x/km/s",
    "v_2y/km/s",
    "v_2z/km/s",
    "v_2/km/s",
    "radius_2/AU",
    "q",
    "projected_separation/AU",
]
data = np.stack([
    df_properties["mass_1/Msun"],
    *binary.primary.state(df_properties["true_anomaly"]).T,
    binary.primary.radius(df_properties["true_anomaly"]),
    binary.primary.speed(df_properties["true_anomaly"]),
    df_properties["mass_1/Msun"]*df_properties["mass_ratio"],
    *binary.secondary.state(df_properties["true_anomaly"]).T,
    binary.secondary.radius(df_properties["true_anomaly"]),
    binary.secondary.speed(df_properties["true_anomaly"]),
    df_properties["mass_ratio"],
    np.sqrt(
        (
            binary.secondary.state(df_properties["true_anomaly"])[:,0]
            - binary.primary.state(df_properties["true_anomaly"])[:,0]
        )**2.
        + (
            binary.secondary.state(df_properties["true_anomaly"])[:,1]
            - binary.primary.state(df_properties["true_anomaly"])[:,1]
        )**2.
    )
])
df_kinematics = pd.DataFrame(data.T, columns=cols)
df_kinematics.to_csv(data_dir + "binary_kinematics.csv", index=False)
# df_kinematics = pd.read_csv(data_dir + "binary_kinematics.csv")

########################################################################
# Compute additional velocity dispersion (exc. spectroscopic binaries)
########################################################################
alpha = np.linspace(0., 1.)
bff = 2.*alpha/(1. + alpha)
sigma2 = 0.5*(
    np.average(
        binary.primary.speed(df_properties["true_anomaly"])**2.,
        weights=weight
    )
    + np.average(
        binary.secondary.speed(df_properties["true_anomaly"])**2.,
        weights=weight
    )
)
delta_sigma2 = bff*sigma2

cols = ["binary_fraction", "additional_dispersion/km/s"]
data = np.array([alpha, delta_sigma2])
df_additional_dispersion = pd.DataFrame(data.T, columns=cols)
df_additional_dispersion.to_csv(
    data_dir + "additional_dispersion.csv", index=False
)

####################################################################
# Visual and spectroscopic contributions
####################################################################
def _sigma2_visual(primary_speed, secondary_speed, weight):
    # Compute visual mean-square inner LOS speed
    if len(primary_speed) == 0:
        res = 0.
    else:
        # res = 0.5*(
        #     np.average(primary_speed**2., weights=weight)
        #     + np.average(secondary_speed**2., weights=weight)
        # )
        try:
            res = 0.5*(
                np.average(primary_speed**2., weights=weight)
                + np.average(secondary_speed**2., weights=weight)
            )
        except ZeroDivisionError:
            # Weights sum to zero
            warnings.warn("weights sum to zero.", UserWarning)
            res = 0.5*(
                np.average(primary_speed**2.)
                + np.average(secondary_speed**2.)
            )

    return res

def _sigma2_spectroscopic(primary_speed, secondary_speed,
                          luminosity_1, luminosity_2,
                          weight):
    # Compute spectroscopic mean-square inner LOS speed
    if len(primary_speed) == 0:
        res = 0.
    else:
        num = luminosity_1*primary_speed + luminosity_2*secondary_speed
        denom = luminosity_1 + luminosity_2
        v_z_weighted = num/denom
        res = np.average(v_z_weighted**2., weights=weight)

    return res

n_bins = 100
s_min = np.logspace(-6., 8., n_bins) # Units: AU
primary_los_speed = binary.primary.state(
    df_properties["true_anomaly"]
)[:,-1]
secondary_los_speed = binary.secondary.state(
    df_properties["true_anomaly"]
)[:,-1]
delta_sigma2 = np.zeros([n_bins, n_bins])
sigma2_visual = np.zeros(n_bins)
sigma2_spectroscopic = np.zeros(n_bins)
for i, s_min_i in enumerate(s_min):
    luminosity_1 = (
        df_properties["luminosity_1/Lsun"].to_numpy()
    )
    luminosity_2 = (
        df_properties["luminosity_2/Lsun"].to_numpy()
    )
    idx_spectroscopic = (
        df_kinematics["projected_separation/AU"] <= s_min_i
    )
    idx_visual = ~idx_spectroscopic
    sigma2_visual[i] = _sigma2_visual(
        primary_los_speed[idx_visual],
        secondary_los_speed[idx_visual],
        weight[idx_visual],
    )
    sigma2_spectroscopic[i] = _sigma2_spectroscopic(
        primary_los_speed[idx_spectroscopic],
        secondary_los_speed[idx_spectroscopic],
        luminosity_1[idx_spectroscopic],
        luminosity_2[idx_spectroscopic],
        weight[idx_spectroscopic],
    )

cols = [
    "spatial_resolution/AU",
    "meansquare_speed/km^2/s^2",
]
data = np.array([s_min, sigma2_visual])
df_visualadditional_dispersion = pd.DataFrame(data.T, columns=cols)
df_visualadditional_dispersion.to_csv(
    data_dir + "visual_meansquare_observed_los_speed.csv",
    index=False
)

cols = [
    "spatial_resolution/AU",
    "meansquare_speed/km^2/s^2",
]
data = np.array([s_min, sigma2_spectroscopic])
df_spectroscopic_dispersion = pd.DataFrame(data.T, columns=cols)
df_spectroscopic_dispersion.to_csv(
    data_dir + "spectroscopic_meansquare_observed_los_speed.csv",
    index=False
)

########################################################################
# Compute additional velocity dispersion (inc. spectroscopic binaries)
########################################################################
alpha = np.linspace(0., 1., n_bins)
s_min = np.logspace(-6., 8., n_bins) # Units: AU

n_sample = 1_000_000
primary_los_speed = binary.primary.state(
    df_properties["true_anomaly"]
)[:n_sample,-1]
secondary_los_speed = binary.secondary.state(
    df_properties["true_anomaly"]
)[:n_sample,-1]
weight = weight[:n_sample]
df_properties = df_properties[:n_sample]
df_kinematics = df_kinematics[:n_sample]

aa = np.zeros([n_bins, n_bins])
ss = np.zeros([n_bins, n_bins])
delta_sigma2 = np.zeros([n_bins, n_bins])
sigma2_visual = np.zeros(n_bins)
sigma2_spectroscopic = np.zeros(n_bins)
for i, alpha_i in enumerate(alpha):
    for j, s_min_j in enumerate(s_min):
        luminosity_1 = (
            df_properties["luminosity_1/Lsun"].to_numpy()
        )
        luminosity_2 = (
            df_properties["luminosity_2/Lsun"].to_numpy()
        )
        idx_spectroscopic = (
            df_kinematics["projected_separation/AU"] <= s_min_j
        )
        idx_visual = ~idx_spectroscopic
        gamma = np.sum(weight[idx_spectroscopic==True])/np.sum(weight)
        # gamma = (
        #     len(idx_spectroscopic[idx_spectroscopic==True])
        #     /len(idx_spectroscopic)
        # )
        sigma2_visual[j] = _sigma2_visual(
            primary_los_speed[idx_visual],
            secondary_los_speed[idx_visual],
            weight[idx_visual],
        )
        sigma2_spectroscopic[j] = _sigma2_spectroscopic(
            primary_los_speed[idx_spectroscopic],
            secondary_los_speed[idx_spectroscopic],
            luminosity_1[idx_spectroscopic],
            luminosity_2[idx_spectroscopic],
            weight[idx_spectroscopic],
        )
        num = (
            2*alpha_i*(1. - gamma)*sigma2_visual[j]
            + alpha_i*gamma*sigma2_spectroscopic[j]
        )
        denom = 1. + alpha_i*(1 - gamma)
        res = num/denom
        if res != res:
            print(alpha_i, s_min_j, idx_spectroscopic)
        aa[i][j] = alpha_i
        ss[i][j] = s_min_j
        delta_sigma2[i][j] = res

cols = [
    "binary_fraction",
    "spatial_resolution/AU",
    "additional_meansquare_speed/km^2/s^2",
]
data = np.array([aa.flatten(), ss.flatten(), delta_sigma2.flatten()])
df_additional_dispersion = pd.DataFrame(data.T, columns=cols)
df_additional_dispersion.to_csv(
    data_dir + "additional_meansquare_observed_los_speed.csv",
    index=False
)
