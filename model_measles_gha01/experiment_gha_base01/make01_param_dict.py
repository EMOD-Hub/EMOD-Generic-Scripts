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

    param_dict[EXP_NAME] = 'Measles-GHA-Base01'
    param_dict[NUM_SIMS] = 600
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

    # Year to end simulation
    P_CON['end_year'] = 2025

    # Coverage of SIAs in WHO calendar
    P_CON['SIA_cover_GHA_2010'] = 0.75
    P_CON['SIA_cover_GHA_2013'] = 0.40
    P_CON['SIA_cover_GHA_2018'] = 0.15

    # Infectivity
    vals = 7.0 + np.random.gamma(30.0, scale=0.133, size=NSIMS)
    P_VAR['R0'] = (np.round(vals, 2)).tolist()

    # R0 seasonality
    P_CON['R0_peak_day'] = 65.0
    P_CON['R0_peak_width'] = 45.0
    P_CON['R0_peak_magnitude'] = 1.3

    # Importation rate per-100k per-day
    P_CON['log10_import_rate'] = -0.50

    # Parameters for gravity model for network connections
    P_CON['net_inf_power'] = 2.0
    P_CON['net_inf_ln_mult'] = 0.2
    P_CON['net_inf_maxfrac'] = 0.1

    # Correlation between acqusition and transmission heterogeneity
    P_CON['corr_acq_trans'] = 0.8

    # Individual risk variance
    P_CON['ind_variance_risk'] = 0.4

    # Base agent weight; less than 10 may have memory issues
    P_CON['agent_rate'] = 25.0

    # Node level overdispersion; 0.0 = Poisson
    P_CON['proc_overdispersion'] = 0.4

    # Reactive campaign case threshold (observed) for admin-1
    P_CON['adm01_case_threshold'] = 1.0e6

    # Reactive campaign case threshold (observed) for admin-1
    P_CON['log10_min_reporting'] = -7

    # Write parameter dictionary
    with open(P_FILE, 'w') as fid01:
        json.dump(param_dict, fid01)

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    write_param_dict()

# *****************************************************************************
