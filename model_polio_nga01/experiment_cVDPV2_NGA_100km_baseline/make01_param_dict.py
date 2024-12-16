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

    param_dict[EXP_NAME] = 'cVDPV2-NGA-100km-base'
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

    #P_VAR['run_number'] = list(range(NSIMS))
    #P_CON['rng_list_offset_yr'] = []
    #P_CON['rng_list_val'] = []

    P_CON['run_number'] = 345
    P_CON['rng_list_offset_yr'] = [3.0, 6.75]
    P_CON['rng_list_val'] = [171, -1]

    # Simulation start / duration
    P_CON['start_year'] = 2018
    P_CON['run_years'] = 6.0

    # Outbreak seeding
    P_CON['seed_location'] = 'AFRO:NIGERIA:JIGAWA:DUTSE'
    P_CON['seed_offset_yr'] = 300.0/365.0

    # Parameters for gravity model for network connections
    P_CON['net_inf_power'] = [2.0]
    P_CON['net_inf_ln_mult'] = [-2.424]

    # Node level overdispersion; 0.0 = Poisson
    P_CON['proc_overdispersion'] = 0.1

    # Correlation between acqusition and transmission heterogeneity
    P_CON['corr_acq_trans'] = 0.8

    # Base agent weight; less than 10 may have memory issues
    P_CON['agent_rate'] = 100.0

    # R0 values for cVDPV, Sabin, nOPV; linear interpolation;
    P_CON['R0'] = 16.0
    P_CON['R0_OPV_mult'] = 0.250
    P_CON['R0_nOPV_mult'] = 0.125

    # Subdivide LGAs into 100km^2 regions
    P_CON['use_10k_res'] = True

    # RI params
    P_CON['ri_start_yr'] = 2030.0

    # Apply the historic SIA calendar; events prior to sim start ignored
    P_CON['sia_calendar'] = True
    P_CON['sia_cutoff'] = 2030.0
    P_CON['sia_coverage'] = 0.5

    # Additional nOPV2 SIAs
    P_CON['nopv2_sia_national'] = []

    # Individual level risk variance (risk of acquisition multiplier;
    # mean = 1.0; log-normal distribution)
    P_CON['ind_variance_risk'] = 3.0

    # Write parameter dictionary
    with open(P_FILE, 'w') as fid01:
        json.dump(param_dict, fid01)

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    write_param_dict()

# *****************************************************************************
