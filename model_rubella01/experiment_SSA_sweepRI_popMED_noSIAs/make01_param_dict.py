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

    param_dict[EXP_NAME] = 'Rubella01-SSA-Demog_Med'
    param_dict[NUM_SIMS] = 3000
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
    P_VAR['R0'] = (1.0 + np.random.gamma(16.0, scale=0.25,
                                         size=NSIMS)).tolist()

    # RI rate for MR
    P_VAR['RI_rate'] = np.random.choice([0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
                                        size=NSIMS).tolist()

    # Context
    P_CON['demog_set'] = 'SSA'
    P_CON['num_nodes'] = 1

    # Log10 of multiplier on exogeneous case importation
    P_CON['log10_import_mult'] = -0.1

    # Use constant vital dynamics
    P_CON['steady_state_demog'] = False

    # Use SIAs
    P_CON['add_campaigns'] = False

    # Initial number of agents
    P_CON['num_agents'] = 100000

    # Write parameter dictionary
    with open(P_FILE, 'w') as fid01:
        json.dump(param_dict, fid01)

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    write_param_dict()

# *****************************************************************************
