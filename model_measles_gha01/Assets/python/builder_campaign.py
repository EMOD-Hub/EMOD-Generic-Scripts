# *****************************************************************************
#
# Campaign file.
#
# *****************************************************************************

import json
import os

import global_data as gdata

import numpy as np

import emod_api.campaign as camp_module

from emod_camp_events import build_node_list, ce_br_force, ce_RI, ce_SIA, \
                             ce_inf_force, ce_inf_mod
from emod_constants import CAMP_FILE

# *****************************************************************************


def campaignBuilder():

    # Variables for this simulation
    PEAK_SIZE = gdata.var_params['R0_peak_magnitude']
    PEAK_TIME = gdata.var_params['R0_peak_day']
    PEAK_WIDE = gdata.var_params['R0_peak_width']

    NODE_DICT = gdata.node_idval
    DAY_MIN = 365.0*(gdata.start_year-gdata.base_year)

    # Note: campaign module itself is the file object; no Campaign class
    ALL_NODES = gdata.demog_object.node_ids

    # Time varying birth rate
    BR_MULT_X = gdata.brate_mult_x
    BR_MULT_Y = gdata.brate_mult_y
    start_day = 365.0*(gdata.start_year-gdata.base_year)
    camp_event = ce_br_force(ALL_NODES, BR_MULT_X, BR_MULT_Y, start_day)
    camp_module.add(camp_event)

    # Add MCV1 RI
    with open(os.path.join('Assets', 'data', 'GHA_MCV1.json')) as fid01:
        mcv1_dict = json.load(fid01)

    time_vec = np.array(mcv1_dict['timevec']) - DAY_MIN
    name_vec = mcv1_dict['namevec']
    mcv1_mat = np.array(mcv1_dict['coverage'])
    for k1 in range(len(name_vec)):
        reg_name = name_vec[k1]
        mcv1_vec = mcv1_mat[k1, :]
        n_list = build_node_list([reg_name], NODE_DICT)

        if (np.amin(time_vec) <= 0.0):
            init_mcv1 = np.interp(0.0, time_vec, mcv1_vec)
        else:
            init_mcv1 = np.mean(mcv1_vec[:3])

        time_list = [0.0] + (time_vec[time_vec > 0.0]).tolist() + \
                    [365.0*100]
        mcv1_list = [init_mcv1] + (mcv1_vec[time_vec > 0.0]).tolist() + \
                    [np.mean(mcv1_vec[-3:])]

        camp_event = ce_RI(n_list, start_day=DAY_MIN,
                           coverage_x=time_list, coverage_y=mcv1_list)
        camp_module.add(camp_event)

    # Add MCV SIAs
    with open(os.path.join('Assets', 'data', 'GHA_MCV_SIA.json')) as fid01:
        dict_sia = json.load(fid01)

    for sia_name in dict_sia:
        sia_day = dict_sia[sia_name]['date']
        if (sia_day < DAY_MIN):
            continue

        SIA_COVER = gdata.var_params['SIA_cover_{:s}'.format(sia_name)]
        age_yr_max = dict_sia[sia_name]['age_yr_max']
        age_yr_min = dict_sia[sia_name]['age_yr_min']

        camp_event = ce_SIA(ALL_NODES, start_day=sia_day, coverage=SIA_COVER,
                            yrs_min=age_yr_min, yrs_max=age_yr_max)
        camp_module.add(camp_event)

    # Add seasonality
    start_day = 365.0*(gdata.start_year-gdata.base_year)
    camp_event = ce_inf_force(ALL_NODES, PEAK_TIME, PEAK_WIDE, PEAK_SIZE,
                              start_day=start_day, dt=gdata.t_step_days)
    camp_module.add(camp_event)

    # Add infectivity trough
    start_day = 365.0*(2020.0-gdata.base_year)
    camp_event = ce_inf_mod(ALL_NODES, start_day=start_day,
                            dt_days=365.0*2.7, mult_val=0.60)
    camp_module.add(camp_event)

    # End file construction
    camp_module.save(filename=CAMP_FILE)

    return None

# *****************************************************************************
