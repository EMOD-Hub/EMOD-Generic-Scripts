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

from emod_camp_events import ce_import_pressure, ce_OPV_SIA, \
                             ce_random_numbers, ce_RI
from emod_constants import CAMP_FILE

# *****************************************************************************


def campaignBuilder():

    # Variables for this simulation
    START_YEAR = gdata.var_params['start_year']
    SIA_CALENDAR = gdata.var_params['sia_calendar']
    SIA_STOP = gdata.var_params['sia_cutoff']
    SIA_COVER = gdata.var_params['sia_coverage']
    SEED_LOCATION = gdata.var_params['seed_location']
    SEED_OFFSET = gdata.var_params['seed_offset_yr']
    RNG_LIST = gdata.var_params['rng_list_offset_yr']
    RNG_VAL = gdata.var_params['rng_list_val']
    RI_START_YR = gdata.var_params['ri_start_yr']
    NODE_DICT = gdata.demog_node

    node_opts = list(NODE_DICT.keys())

    # Note: campaign module itself is the file object; no Campaign class
    ALL_NODES = gdata.demog_object.node_ids

    # Use SIA calendar for OPV2 schedule
    if (SIA_CALENDAR):
        with open(os.path.join('Assets', 'data', 'sia_dat_NGA.json')) as fid01:
            dict_sia = json.load(fid01)

        for sia_name in dict_sia:
            sia_obj = dict_sia[sia_name]

            startday = sia_obj['date']
            if (startday > SIA_STOP):
                continue

            if (sia_obj['type'] == 'sabin2'):
                clade = 0
                genome = gdata.boxes_nopv2
            elif (sia_obj['type'] == 'nopv2'):
                clade = 1
                genome = 0

            node_list = list()
            for targ_val in sia_obj['nodes']:
                for nname in node_opts:
                    if ((nname == targ_val) or
                       (nname.startswith(targ_val+':'))):
                        node_list.append(NODE_DICT[nname])

            camp_event = ce_OPV_SIA(node_list, start_day=startday,
                                    coverage=SIA_COVER,
                                    clade=clade, genome=genome)
            camp_module.add(camp_event)

    # Seed infections
    node_list = list()
    for nname in node_opts:
        if ((nname == SEED_LOCATION) or (nname.startswith(SEED_LOCATION+':'))):
            node_list.append(NODE_DICT[nname])

    # Preserve size of outbreak; select single node for initial location
    node_list = [node_list[-1]]
    cvdpv_gen = gdata.boxes_nopv2+gdata.boxes_sabin2
    startday = 365.0*(START_YEAR-gdata.base_year+SEED_OFFSET)
    camp_event = ce_import_pressure(node_list, start_day=startday,
                                    genome=cvdpv_gen,
                                    duration=gdata.seed_inf_dt,
                                    magnitude=gdata.seed_inf_num)
    camp_module.add(camp_event)

    # Add SIAs
    #for syear in [7.5, 8.5]:
    #    start_day = 365.0*(START_YEAR-gdata.base_year+syear)
    #    camp_event = ce_OPV_SIA(ALL_NODES, start_day=start_day,
    #                            coverage=SIA_COVER, clade=1, genome=0)
    #    camp_module.add(camp_event)

    ## Add RI
    start_day = 365.0*(START_YEAR-gdata.base_year+6.5)
    camp_event = ce_RI(ALL_NODES, [0.0], [0.5], start_day=start_day,
                       base_take=0.7, age_one=180.0, age_std=5.0)
    camp_module.add(camp_event)


    # Random number stream offset
    for (yr_off, nval) in zip(RNG_LIST, RNG_VAL):
        startday = 365.0*(START_YEAR-gdata.base_year+yr_off)
        if (nval < 0):
            nval = gdata.sim_index
        camp_event = ce_random_numbers(ALL_NODES, start_day=startday,
                                       numbers=nval)
        camp_module.add(camp_event)

    # End file construction
    camp_module.save(filename=CAMP_FILE)

    return None

# *****************************************************************************
