#!/usr/bin/env python

"""Plot binary properties"""

import pathlib
import warnings
import yaml
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import plot

mpl.style.use("sm")

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

population = config["population"]
n_bins = 100

data_dir = population + "/data/"
df_properties = pd.read_csv(data_dir + "binary_properties.csv")
try:
    weight = df_properties["weight"].values
except KeyError:
    warnings.warn(
        "data frame does not contain a weight column.", UserWarning
    )
    df_properties["weight"] = 1.
    weight = df_properties["weight"].values

figures_dir = population + "/figures/"
try:
    pathlib.Path(figures_dir).mkdir(parents=True)
except FileExistsError:
    pass

fig0, ax0 = plot.plot()
ax0.hist(
    df_properties["mass_1/Msun"], bins=np.logspace(-2., 2., n_bins),
    histtype="step", density=True, weights=weight
)
ax0.set_ylim(1.e-4, 1.e1)
ax0.set_xscale("log")
ax0.set_yscale("log")
ax0.set_xlabel(r"$M_{1}/\mathrm{M}_{\odot}$")
ax0.set_ylabel(r"$\hat{f}$")
fig0.savefig(figures_dir + "01_primary_mass.jpg")
fig0.savefig(figures_dir + "01_primary_mass.pdf")

fig2, ax2 = plot.plot()
ax2.hist(
    df_properties["mass_1/Msun"], bins=np.logspace(-2., 2., n_bins),
    histtype="step", density=True, weights=weight
)
ax2.set_ylim(1.e-4, 1.e1)
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel(r"$M_{2}/\mathrm{M}_{\odot}$")
ax2.set_ylabel(r"$\hat{f}$")
fig2.savefig(figures_dir + "02_secondary_mass.jpg")
fig2.savefig(figures_dir + "02_secondary.pdf")

fig3, ax3 = plot.plot()
ax3.hist(
    df_properties[["mass_1/Msun", "mass_2/Msun"]].values.ravel(),
    bins=np.logspace(-2., 2., n_bins), histtype="step", density=True,
    weights=np.hstack([weight, weight])
)
ax3.set_ylim(1.e-4, 1.e1)
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlabel(r"$M/\mathrm{M}_{\odot}$")
ax3.set_ylabel(r"$\hat{f}$")
fig3.savefig(figures_dir + "03_all_mass.jpg")
fig3.savefig(figures_dir + "03_al_mass.pdf")

fig4, ax4 = plot.plot()
ax4.hist(
    df_properties["mass_ratio"], bins=np.linspace(-0.2, 1.2, n_bins),
    histtype="step", density=True, weights=weight
)
ax4.set_xlabel(r"$q$")
ax4.set_ylabel(r"$\hat{f}$")
fig4.savefig(figures_dir + "04_mass_ratio.jpg")
fig4.savefig(figures_dir + "04_mass_ratio.pdf")

fig5, ax5 = plot.plot()
ax5.hist(
    np.log10(df_properties["semimajor_axis/AU"]),
    bins=np.linspace(-6., 8., n_bins), histtype="step", density=True,
    weights=weight
)
ax5.set_xlabel(r"$a/\mathrm{AU}$")
ax5.set_ylabel(r"$\hat{f}$")
fig5.savefig(figures_dir + "05_semimajor_axis.jpg")
fig5.savefig(figures_dir + "05_semimajor_axis.pdf")

fig6, ax6 = plot.plot()
ax6.hist(
    np.log10(df_properties["period/day"]), bins=np.linspace(-3., 13., n_bins),
    histtype="step", density=True, weights=weight
)
ax6.set_xlabel(r"$\log(P/\mathrm{day})$")
ax6.set_ylabel(r"$\hat{f}$")
fig6.savefig(figures_dir + "06_period.jpg")
fig6.savefig(figures_dir + "06_period.pdf")

fig7, ax7 = plot.plot()
ax7.hist(
    df_properties["eccentricity"], bins=np.linspace(-0.2, 1.2, n_bins),
    histtype="step", density=True, weights=weight
)
ax7.set_xlim(-0.2, 1.2)
ax7.set_xlabel(r"$e$")
ax7.set_ylabel(r"$\hat{f}$")
fig7.savefig(figures_dir + "07_eccentricity.jpg")
fig7.savefig(figures_dir + "07_eccentricity.pdf")

fig8, ax8 = plot.plot()
ax8.hist(
    df_properties["true_anomaly"], bins=n_bins, histtype="step", density=True,
    weights=weight
)
ax8.set_xlim(-np.pi/2., 5.*np.pi/2.)
ax8.set_xlabel(r"$\theta$")
ax8.set_ylabel(r"$\hat{f}$")
fig8.savefig(figures_dir + "08_true_anomaly.jpg")
fig8.savefig(figures_dir + "08_true_anomaly.pdf")

fig9, ax9 = plot.plot()
ax9.hist(
    df_properties["longitude_of_ascending_node"], bins=n_bins, histtype="step",
    density=True, weights=weight
)
ax9.set_xlim(-np.pi/2., 5.*np.pi/2.)
ax9.set_xlabel(r"$\Omega$")
ax9.set_ylabel(r"$\hat{f}$")
fig9.savefig(figures_dir + "09_longitude_of_ascending_node.jpg")
fig9.savefig(figures_dir + "09_longitude_of_ascending_node.pdf")

fig10, ax10 = plot.plot()
ax10.hist(
    df_properties["inclination"], bins=n_bins, histtype="step", density=True,
    weights=weight
)
ax10.set_xlim(-np.pi/2., 3.*np.pi/2.)
ax10.set_xlabel(r"$i$")
ax10.set_ylabel(r"$\hat{f}$")
fig10.savefig(figures_dir + "10_inclination.jpg")
fig10.savefig(figures_dir + "10_inclination.pdf")

fig11, ax11 = plot.plot()
ax11.hist(
    df_properties["argument_of_pericentre"], bins=n_bins, histtype="step",
    density=True, weights=weight
)
ax11.set_xlim(-np.pi/2., 5.*np.pi/2.)
ax11.set_xlabel(r"$\omega$")
ax11.set_ylabel(r"$\hat{f}$")
fig11.savefig(figures_dir + "11_argument_of_pericentre.jpg")
fig11.savefig(figures_dir + "11_argument_of_pericentre.pdf")

plt.close("all")
