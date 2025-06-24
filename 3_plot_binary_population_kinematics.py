#!/usr/bin/env python

"""Plot binary kinematics"""

import pathlib
import warnings
import yaml
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import dyad
import plot

mpl.style.use("sm")

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

population = config["population"]
n_bins = 100

data_dir = population + "/data/"
df_kinematics = pd.read_csv(data_dir + "binary_kinematics.csv")
df_properties = pd.read_csv(data_dir + "binary_properties.csv")
try:
    weight = df_properties["weight"]
except KeyError:
    warnings.warn(
        "data frame does not contain a weight column.", UserWarning
    )
    df_properties["weight"] = 1.
    weight = df_properties["weight"]

figures_dir = population + "/figures/"
try:
    pathlib.Path(figures_dir).mkdir(parents=True)
except FileExistsError:
    pass
    
# Speeds of binary components
bins = np.linspace(-4., 4., n_bins)
hist_0 = np.histogram(
    np.log10(df_kinematics["v_1/km/s"]), bins=bins, density=True,
    weights=weight
)
hist_1 = np.histogram(
    np.log10(df_kinematics["v_2/km/s"]), bins=bins, density=True,
    weights=weight
)

fig, ax = plot.plot()
ax.stairs(*hist_0, label=r"$i = \mathrm{b}_{1}$", zorder=1_000)
ax.stairs(*hist_1, label=r"$i = \mathrm{b}_{2}$", color="red")
ax.legend(frameon=False)
ax.set_ylim(0., 0.6)
ax.set_xlabel(r"$\log_{10}(v_{i}/\mathrm{km}~\mathrm{s}^{-1})$")
ax.set_ylabel(r"$\hat{f}$")
fig.savefig(figures_dir + "binary_speed.jpg")
fig.savefig(figures_dir + "binary_speed.pdf")

# Projected separation of binary components
galactic_distance = 1.e6 # Unit: pc
s_resolution = [
    galactic_distance*1.e-3,
    galactic_distance*1.e-2,
    galactic_distance*1.e-1
] # Unit: arcsec
s_min = 1.e-6
s_max = 1.e8
bins = np.logspace(-6., 8., n_bins)
# Find counts using unnormalized histogram
hist, edges = np.histogram(
    df_kinematics["projected_separation/AU"], bins=bins,
    weights=weight
)
hist = np.cumsum(hist)
hist = hist/hist[-1]

fig, ax = plot.plot()
ax.plot(edges[:-1], hist)
ax.vlines(s_resolution, 0., 1., ls="dashed")
ax.text(
    s_resolution[0] + 100.,
    0. + 0.05,
    r"$\delta\theta = 0.001~\text{as}$",
    horizontalalignment="left",
    verticalalignment="bottom",
    rotation=90.,
)
ax.text(
    s_resolution[1] + 100.,
    0. + 0.05,
    r"$\delta\theta = 0.01~\text{as}$",
    horizontalalignment="left",
    verticalalignment="bottom",
    rotation=90.,
)
ax.text(
    s_resolution[2] - 100.,
    0. + 0.05,
    r"$\delta\theta = 0.1~\text{as}$",
    horizontalalignment="left",
    verticalalignment="bottom",
    rotation=90.,
)
ax.set_xlim(1.e-6, 1.e8)
ax.set_ylim(0., 1.)
ax.set_xscale("log")
ax.set_xlabel(r"$\xi/\text{au}$")
ax.set_ylabel(r"$\hat{F}_{S}$")
fig.savefig(figures_dir + "on-sky_separation.jpg")
fig.savefig(figures_dir + "on-sky_separation.pdf")

# Additional dispersion
df = pd.read_csv(data_dir + "additional_dispersion.csv")
fig, ax = plot.plot()
ax.plot(df["binary_fraction"], df["additional_dispersion/km/s"])
# ax.set_ylim(0., 1.5e3)
ax.set_xlabel(r"$\alpha$")
ax.set_ylabel(r"$\delta{}\sigma^{2}_{V}/\mathrm{km}^{2}~\mathrm{s}^{-2}$")
fig.savefig(figures_dir + "additional_mean_square_speed.jpg")
fig.savefig(figures_dir + "additional_mean_square_speed.pdf")

# Internal LOS velocity: resolved and unresolved
v_z = np.concatenate(
    [df_kinematics["v_1z/km/s"], df_kinematics["v_2z/km/s"]]
)
bins = np.linspace(-50., 50., n_bins)
hist_0 = np.histogram(
    v_z, bins=bins, density=True#, weights=np.concatenate([weight, weight])
)

v_z = np.average(
    df_kinematics[["v_1z/km/s", "v_2z/km/s"]],
    # weights=df_properties[
    #     ["primary_luminosity/Lsun", "secondary_luminosity/Lsun"]
    # ],
    axis=1
) # Luminosity-weighted average LOS velocity
hist_1 = np.histogram(
    v_z, bins=bins, density=True, weights=weight
)

fig, ax = plot.plot()
ax.stairs(*hist_0, label=r"visual", zorder=1_000)
ax.stairs(*hist_1, label=r"spectroscopic", color="green")
ax.legend(frameon=False)
# ax.set_xlim(-5., 5.)
ax.set_xlim(-25., 25.)
# ax.set_ylim(0., 1.4)
ax.set_ylim(1.e-4, 1.e0)
ax.set_yscale("log")
ax.set_xlabel(r"$v_{z}/\mathrm{km}~\mathrm{s}^{-1}$")
ax.set_ylabel(r"$\hat{f}_{V_{z}}$")
fig.savefig(figures_dir + "resolved_unresolved_los_velocity.jpg")
fig.savefig(figures_dir + "resolved_unresolved_los_velocity.pdf")

# Additional mean-square speed
df = pd.read_csv(
    data_dir + "additional_meansquare_observed_los_speed.csv"
)
aa = df["binary_fraction"].to_numpy().reshape([n_bins, n_bins])
ss = df["spatial_resolution/AU"].to_numpy().reshape([n_bins, n_bins])
delta_sigma2 = df["additional_meansquare_speed/km^2/s^2"]
delta_sigma2 = delta_sigma2.to_numpy().reshape([n_bins, n_bins])

if population == "moe2017_canonical":
    manual = [
        [0.1, 6.],
        [0.35, 6.],
        [0.8, 6.],
        [0.1, -4.5],
        [0.15, -4.5],
        [0.25, -4.5],
        [0.35, -4.5],
        [0.5, -4.5],
        [0.75, -4.5],
    ]
    major_axis_min = 2.*dyad.semimajor_axis_from_period(10.**0.2, 0.8, 0.8)
    levels = (np.arange(0., 50., 2)**2.)
elif population == "moe2017_bottom_light":
    manual = [
        [0.1, 6.],
        [0.6, 6.],
        # [0.05, -4.5],
        [0.1, -4.5],
        [0.15, -4.5],
        [0.25, -4.5],
        [0.4, -4.5],
        [0.6, -4.5],
        [1., -4.5],
    ]
    major_axis_min = 2.*dyad.semimajor_axis_from_period(10.**0.2, 0.8, 0.8)
    levels = (np.arange(0., 50., 2)**2.)[2::2]
elif population == "duquennoy1991":
    manual = [
        [0.1, 6.],
        [0.35, 6.],
        [0.8, 6.],
        [0.1, -4.5],
        [0.15, -4.5],
        [0.25, -4.5],
        [0.35, -4.5],
        [0.5, -4.5],
        [0.75, -4.5],
    ]
    major_axis_min = 2.*dyad.semimajor_axis_from_period(10.**(-2.4), 0.9, 0.9)
    levels = (np.arange(0., 50., 2)**2.)

fig, ax = plot.plot()
im = ax.pcolormesh(
    aa, np.log10(ss), delta_sigma2, rasterized=True, vmin=0.,
    vmax=1.2*np.max(delta_sigma2),
    cmap="Greys"
)
cs = ax.contour(
    aa, np.log10(ss), delta_sigma2, colors="k", levels=levels
)
cb = ax.clabel(
    cs, inline=True, use_clabeltext=True, inline_spacing=5.,
    manual=manual
)
ax.set_xlim(0., 1.)
ax.set_ylim(-6, 8.)
ax.set_xlabel(r"$\alpha$")
ax.set_ylabel(r"$\log_{10}(\xi/\text{au})$")
ax.set_title(r"$\delta\sigma^{2}_{V}/\text{km}^{2}~\text{s}^{-2}$")
ax.hlines(np.log10(major_axis_min), 0., 1., ls="dashed")
fig.savefig(
    figures_dir + "additional_meansquare_observed_los_speed.jpg"
)
fig.savefig(
    figures_dir + "additional_meansquare_observed_los_speed.pdf"
)

# Visual and spectroscopic binary mean-square speed
df_visual = pd.read_csv(
    data_dir + "visual_meansquare_observed_los_speed.csv"
)
df_spectroscopic = pd.read_csv(
    data_dir + "spectroscopic_meansquare_observed_los_speed.csv"
)

s = df_visual["spatial_resolution/AU"]
sigma2_visual = df_visual["meansquare_speed/km^2/s^2"]
sigma2_spectroscopic = df_spectroscopic["meansquare_speed/km^2/s^2"]

fig, ax = plot.plot()
levels = (np.arange(0., 50., 2)**2.)[1:]
ax.plot(np.log10(s), sigma2_visual, label=r"$i = \text{v}$")
ax.plot(np.log10(s), sigma2_spectroscopic, label=r"$i = \text{s}$")
ax.set_ylim(10., 10_000.)
ax.set_yscale("log")
ax.set_xlabel(r"$\log_{10}(\xi/\text{au})$")
ax.set_ylabel(r"$\sigma^{2}_{V_{i, z}}$")
ax.legend(frameon=False)
fig.savefig(
    figures_dir + "visual_and_spectroscopic_meansquare_observed_los_speed.jpg"
)
fig.savefig(
    figures_dir + "visual_and_spectroscopic_meansquare_observed_los_speed.pdf"
)

plt.close("all")
