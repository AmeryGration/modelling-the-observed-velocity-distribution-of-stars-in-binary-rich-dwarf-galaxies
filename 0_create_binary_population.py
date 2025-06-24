#!/usr/bin/env python

"""Compute binary population"""

import time
import yaml
import duquennoy1991
import moe2017_bottom_light
import moe2017_canonical

########################################################################
# Load variables from configuration file
########################################################################
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

population = config["population"]
n_binary = config["n_binary"]

########################################################################
# Create binary population data
########################################################################
create_binary_population = getattr(
    globals()[population], "create_binary_population"
)
print(time.ctime())
create_binary_population(n_binary)
