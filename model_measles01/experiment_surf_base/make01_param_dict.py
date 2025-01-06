# *****************************************************************************
#
# *****************************************************************************

import json
import os
import sys

import numpy as np

# Ought to go in emodpy
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'local_python')))
from py_assets_common.emod_constants import EXP_C, EXP_V, EXP_NAME, NUM_SIMS, \
                                            P_FILE

# *****************************************************************************

# This script makes a json dictionary that is used by the pre-processing script
# in EMOD. Variable names defined here will be available to use in creating
# the input files.

# The pre-process script will open the json dict created by this method. For
# everything with the EXP_V key, that script will assume a list and get a value
# from that list based on the sim index. For everything with the EXP_C key, it
# will assume a single value and copy that value.


def write_param_dict():

    # Setup
    param_dict = dict()

    param_dict[EXP_NAME] = 'Measles01-SurfaceBase'
    param_dict[NUM_SIMS] = 9000
    param_dict[EXP_V] = dict()
    param_dict[EXP_C] = dict()

    # Random number consistency
    np.random.seed(4)

    # Convenience naming
    NSIMS = param_dict[NUM_SIMS]
    P_VAR = param_dict[EXP_V]
    P_CON = param_dict[EXP_C]

    # Run number (EMOD random seed)
    P_VAR['run_number'] = list(range(NSIMS))

    # Infectivity
    vals = 10.0 + np.random.gamma(30.0, scale=0.133, size=NSIMS)
    P_VAR['R0'] = (np.round(vals, 2)).tolist()

    # RI params
    vals = np.random.choice(np.arange(0.2, 1.01, 0.04), size=NSIMS)
    P_VAR['MCV1'] = (np.round(vals, 2)).tolist()
    vals = 365.0*np.random.choice(np.arange(3, 19)/12, size=NSIMS)
    P_VAR['MCV1_age'] = (np.round(vals, 2)).tolist()

    # Reference year for population; uses UN WPP for DRC
    P_CON['ref_year'] = 2040

    # SIA parameters
    P_CON['sia_start_year'] = 1000
    P_CON['sia_coverage'] = 0.60
    P_CON['sia_min_age_yr'] = 0.75

    # Log10 of multiplier on exogeneous case importation
    P_CON['log10_import_mult'] = 1.0

    # RI params
    P_CON['Age_Take'] = True
    P_CON['MCV2'] = 0.0

    # Maternal protection effectiveness
    P_CON['mat_factor'] = 1.0

    # Write parameter dictionary
    with open(P_FILE, 'w') as fid01:
        json.dump(param_dict, fid01)

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    write_param_dict()

# *****************************************************************************
