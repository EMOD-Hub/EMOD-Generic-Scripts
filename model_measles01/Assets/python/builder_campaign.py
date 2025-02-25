# *****************************************************************************
#
# Campaign file.
#
# *****************************************************************************

import global_data as gdata

import numpy as np

import emod_api.campaign as camp_module

from emod_camp_events import ce_br_force, ce_inf_force, ce_RI, ce_SIA
from emod_constants import CAMP_FILE

# *****************************************************************************


def campaignBuilder():

    # Variables for this simulation
    MCV1_RATE = gdata.var_params['MCV1']
    MCV2_FRAC = gdata.var_params['MCV2']

    MCV1_AGE = gdata.var_params['MCV1_age']
    AGE_DEP = gdata.var_params['Age_Take']

    SIA_START = gdata.var_params['sia_start_year']
    SIA_COVERAGE = gdata.var_params['sia_coverage']
    SIA_MIN_AGE_YR = gdata.var_params['sia_min_age_yr']

    MAT_FACTOR = gdata.var_params['mat_factor']

    # Note: campaign module itself is the file object; no Campaign class
    ALL_NODES = gdata.demog_object.node_ids

    # Time varying birth rate
    BR_MULT_X = gdata.brate_mult_x
    BR_MULT_Y = gdata.brate_mult_y
    start_day = 365.0*(gdata.start_year-gdata.base_year)
    camp_event = ce_br_force(ALL_NODES, BR_MULT_X, BR_MULT_Y, start_day)
    camp_module.add(camp_event)

    # Add seasonality
    start_day = 365.0*(gdata.start_year-gdata.base_year)
    camp_event = ce_inf_force(ALL_NODES, 15.0, 60.0, 1.30,
                              dt=gdata.t_step_days)
    camp_module.add(camp_event)

    # Add MCV
    start_day = 365.0*(gdata.start_year-gdata.base_year)
    acq_fact = MAT_FACTOR/2.0
    camp_event = ce_RI(ALL_NODES, coverage=MCV1_RATE, start_day=start_day,
                       base_take=0.95, acq_fact=acq_fact, age_dep=AGE_DEP,
                       age_one=MCV1_AGE, frac_two=MCV2_FRAC, age_std=60.0)
    camp_module.add(camp_event)

    # Add SIAs
    start_day = 365.0*(gdata.start_year-gdata.base_year)
    sia_year = SIA_START
    sia_rate = 1.0/(1.0-MCV1_RATE+0.001)
    acq_fact = MAT_FACTOR/2.0

    while (sia_year < gdata.run_years):
        sia_year = sia_year + max(2.0, np.random.poisson(sia_rate))
        start_sia = 365.0*sia_year+start_day

        camp_event = ce_SIA(ALL_NODES, start_day=start_sia,
                            yrs_min=SIA_MIN_AGE_YR, coverage=SIA_COVERAGE,
                            base_take=0.95, acq_fact=acq_fact, age_dep=AGE_DEP)
        camp_module.add(camp_event)

    # End file construction
    camp_module.save(filename=CAMP_FILE)

    return None

# *****************************************************************************
