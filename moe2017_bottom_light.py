#!/usr/bin/env python

"""Create binary population of ZAMS stars with bottom-light IMF"""

import sys
import pathlib
import time
import numpy as np
import scipy as sp
import pandas as pd

import dyad
import dyad.stats as stats
import zams_hrd

def create_binary_population(n_binary=1_000_000):
    # Sample primary mass
    rv_mass = sp.stats.truncpareto(0.3, (40. - 0.)/0.8, scale=0.8)
    masses_primary = rv_mass.rvs(size=n_binary)

    # Sample period
    rv_logperiod = stats.log_period.moe2017(masses_primary)
    log10_periods = rv_logperiod.rvs()
    periods = 10.**log10_periods

    # Sample mass ratio
    rv_mass_ratio = stats.mass_ratio.moe2017(
        log10_periods, masses_primary
    )
    mass_ratios = rv_mass_ratio.rvs()

    # Computer secondary masses
    masses_secondary = masses_primary*mass_ratios

    # Convert period to semimajor axis
    semimajor_axes = dyad.semimajor_axis_from_period(
        periods, masses_primary, masses_primary*mass_ratios,
    )

    # Sample eccentricity
    mask = np.isfinite(
        stats.eccentricity.moe2017(log10_periods, masses_primary).support()[0]
    )
    rv_eccentricity = stats.eccentricity.moe2017(
        log10_periods[mask], masses_primary[mask]
    )
    eccentricities = np.zeros(n_binary)
    eccentricities[mask] = rv_eccentricity.rvs()

    # Sample true anomaly
    rv_true_anomaly = stats.true_anomaly(eccentricities)
    true_anomalies = rv_true_anomaly.rvs()

    # Sample longitude of the ascending node
    rv_longitude_of_ascending_node = (
        stats.longitude_of_ascending_node()
    )
    longitudes_of_ascending_node = (
        rv_longitude_of_ascending_node.rvs(size=n_binary)
    )

    # Sample longitude of the inclination
    rv_inclination = stats.inclination()
    inclinations = rv_inclination.rvs(size=n_binary)

    # Sample longitude of the argument of periapsis
    rv_argument_of_pericentre = stats.argument_of_pericentre()
    arguments_of_pericentre = rv_argument_of_pericentre.rvs(size=n_binary)

    # Luminosity
    primary_luminosities = zams_hrd.luminosity(masses_primary)
    secondary_luminosities = zams_hrd.luminosity(masses_secondary)

    # Effective temperature
    primary_temperatures = zams_hrd.effective_temperature(masses_primary)
    secondary_temperatures = zams_hrd.effective_temperature(masses_secondary)

    # Radius
    primary_radii = zams_hrd.radius(masses_primary)
    secondary_radii = zams_hrd.radius(masses_secondary)

    # Save data in CSV form
    cols = [
        "mass_1/Msun",
        "mass_2/Msun",
        "mass_ratio",
        "semimajor_axis/AU",
        "true_anomaly",
        "eccentricity",
        "longitude_of_ascending_node",
        "inclination",
        "argument_of_pericentre",
        "period/day",
        "luminosity_1/Lsun",
        "luminosity_2/Lsun",
        "effective_temperature_1/K",
        "effective_temperature_2/K",
        "radius_1/Rsun",
        "radius_2/Rsun",
    ]
    data = np.array([
        masses_primary,
        masses_secondary,
        mass_ratios,
        semimajor_axes,
        true_anomalies,
        eccentricities,
        longitudes_of_ascending_node,
        inclinations,
        arguments_of_pericentre,
        periods,
        primary_luminosities,
        secondary_luminosities,
        primary_temperatures,
        secondary_temperatures,
        primary_radii,
        secondary_radii,
    ])
    df = pd.DataFrame(data.T, columns=cols)

    data_dir = "moe2017_bottom_light/data/"
    try:
        pathlib.Path(data_dir).mkdir(parents=True)
    except FileExistsError:
        pass
    df.to_csv(data_dir + "5_binary_properties.csv", index=False)
