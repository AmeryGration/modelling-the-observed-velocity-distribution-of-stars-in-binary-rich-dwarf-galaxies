# Modelling the observed velocity distribution of stars in binary-rich dwarf galaxies

This repository contains the code used in the analysis of binary-rich galaxies performed by Gration et al. [[1]](#1). That analysis quantified the difference in the line-of-sight (LOS) velocity dispersions for binary-free and binary-rich galaxies (called the `additional mean-square LOS velocity') under the assumption that the binary-star population consisted of either
- zero-age main-sequence stars with primary-star masses distributed according to the canonical initial mass function,
- zero-age main-sequence stars with primary-star masses distributed according to a bottom-light initial mass function, or
- present-day solar-type stars.
This analysis was done by realizing a sample of a binary-star population using Monte Carlo methods.

## Installation

Requirements file.

## Usage

To run the analysis first modify the configuration file, `config.yaml`, in order to specify the population type, `population` (one of `duquennoy1991`, `moe2017_canonical`, or `moe2017_bottom_light`), and population size, `n_binary` (an integer).

Then run the following scripts in order:
- `0_create_binary_population.py`, which synthesizes a population of binary stars,
- `1_plot_binary_population_properties.py`, which plots histograms of the physical properties of these binary stars,
- `2_compute_binary_population_kinematics.py`, which computes the kinematic properties of each component of each binary star along with the additional LOS velocity dispersion as a function of binary fraction and spatial resolution, and
- `3_plot_binary_population_kinematics.py`, which plots histograms of the binary star population's kinematic properties and the additional LOS velocity dispersion. 

Executing the script `0_create_binary_population.py` in turn executes one of the scripts `duquennoy1991.py`, `moe2017_canonical.py`, or `moe2017_bottom_light.py` according to the value of the key `population` in `config.yaml`. All data are saved to `population/data/`. All plots are saved to `population/figures/`.

## References
<a id="1">[1]</a> 
Gration, A., Hendriks, D. D., Das, P., Heber, D., and R.G. Izzard (2025). `Modelling the observed velocity distribution of stars in binary-rich ultra-faint dwarf galaxies' in _Monthly Notices of the Royal Astronomical
Society_.
