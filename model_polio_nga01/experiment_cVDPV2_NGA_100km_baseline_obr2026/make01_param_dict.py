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

    param_dict[EXP_NAME] = 'cVDPV2-NGA-100km-base-y18-obr2026'
    param_dict[NUM_SIMS] = 360
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

    P_CON['run_number'] = 2611
    P_CON['rng_list_offset_yr'] = [1.2, 1.3, 1.7, 2.5, 3.6, 4.1, 6.0, 7.5, 8.4]
    P_CON['rng_list_val'] = [223, 14, 337, 170, 158, 652, 32, 297, -1]

    # Simulation start / duration
    P_CON['start_year'] = 2017
    P_CON['run_years'] = 18.0

    # Parameters for gravity model for network connections
    P_CON['net_inf_power'] = [1.5]
    P_CON['net_inf_ln_mult'] = [-2.3]

    # Node level overdispersion; 0.0 = Poisson
    P_CON['proc_overdispersion'] = 0.4

    # Base agent weight; less than 10 may have memory issues
    P_CON['agent_rate'] = 20.0

    # R0 values for cVDPV, Sabin, nOPV; linear interpolation;
    P_CON['R0'] = 18.0
    P_CON['R0_OPV_mult'] = 0.250
    P_CON['R0_nOPV_mult'] = 0.125
    P_CON['R0_sig_scale'] = 32.0
    P_CON['R0_min_mult'] = 0.2
    P_CON['R0_mods'] = [
                       ('AFRO:NIGERIA:BORNO', 2017.5, 20.0, 1.5),
                       ('AFRO:CHAD', 2020.0, 1.0, 1.5),
                       ('AFRO:NIGERIA:KANO', 2017.5, 3.0, 0.4),
                       ('AFRO:NIGERIA:KANO', 2020.5, 0.5, 0.4),
                       ('AFRO:NIGERIA:KANO', 2021.0, 0.4, 0.4),
                       ('AFRO:NIGER', 2021.3, 0.5, 1.5),
                       ('AFRO:NIGERIA:KEBBI', 2021.3, 0.5, 1.5),
                       ('AFRO:NIGERIA:SOKOTO', 2021.3, 0.5, 1.5),
                       ]

    # Individual level risk variance (risk of acquisition multiplier;
    # mean = 1.0; log-normal distribution)
    P_CON['ind_variance_risk'] = 4.02

    # Subdivide LGAs into 100km^2 regions
    P_CON['use_10k_res'] = True

    # RI params
    P_CON['ri_start_yr'] = 2100

    # Apply SIA calendars
    P_CON['sia_end_yr'] = 2100
    P_CON['sia_base_coverage'] = 0.42
    P_CON['sia_coverage_scale'] = 1.2
    P_CON['sia_base_vax_take'] = 0.75
    P_CON['sia_plan_file'] = ''

    # Outbreak response SIAs
    P_CON['use_obr'] = True
    P_CON['obr_start'] = (2026-1900)*365.0

    # Write parameter dictionary
    with open(P_FILE, 'w') as fid01:
        json.dump(param_dict, fid01)

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    write_param_dict()

# *****************************************************************************
