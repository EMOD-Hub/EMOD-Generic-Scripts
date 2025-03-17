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

    param_dict[EXP_NAME] = 'cVDPV2-NGA-100km-base-RI-SIAs'
    param_dict[NUM_SIMS] = 400
    param_dict[EXP_V] = dict()
    param_dict[EXP_C] = dict()

    # Random number consistency
    np.random.seed(4)

    # Convenience naming
    NSIMS = param_dict[NUM_SIMS]
    P_VAR = param_dict[EXP_V]
    P_CON = param_dict[EXP_C]

    # Run number (EMOD random seed)
    #P_VAR['run_number'] = list(range(NSIMS))
    #P_CON['rng_list_offset_yr'] = []
    #P_CON['rng_list_val'] = []

    P_CON['run_number'] = 526
    P_CON['rng_list_offset_yr'] = [1.3, 2.0, 2.9, 4.5, 7.0, 7.8]
    P_CON['rng_list_val'] = [451, 519, 418, 419, 39, -1]

    # Simulation start / duration
    P_CON['start_year'] = 2017
    P_CON['run_years'] = 18.0

    # Parameters for gravity model for network connections
    P_CON['net_inf_power'] = [2.0]
    P_CON['net_inf_ln_mult'] = [-2.424]

    # Node level overdispersion; 0.0 = Poisson
    P_CON['proc_overdispersion'] = 0.4

    # Correlation between acqusition and transmission heterogeneity
    P_CON['corr_acq_trans'] = 0.8

    # Base agent weight; less than 10 may have memory issues
    P_CON['agent_rate'] = 25.0

    # R0 values for cVDPV, Sabin, nOPV; linear interpolation;
    P_CON['R0'] = 14.0
    P_CON['R0_OPV_mult'] = 0.250
    P_CON['R0_nOPV_mult'] = 0.125
    P_CON['R0_sig_scale'] = 24.0
    P_CON['R0_min_mult'] = 0.2

    # Individual level risk variance (risk of acquisition multiplier;
    # mean = 1.0; log-normal distribution)
    P_CON['ind_variance_risk'] = 4.0

    # Subdivide LGAs into 100km^2 regions
    P_CON['use_10k_res'] = True

    # RI params
    P_CON['ri_start_yr'] = 2025.0

    # Apply SIA calendars
    P_CON['sia_base_coverage'] = 0.4
    P_CON['sia_coverage_scale'] = 0.8
    P_CON['sia_base_vax_take'] = 0.7
    P_CON['sia_plan_file'] = 'sia_plan_ALL_2Y.csv'

    # Write parameter dictionary
    with open(P_FILE, 'w') as fid01:
        json.dump(param_dict, fid01)

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    write_param_dict()

# *****************************************************************************
