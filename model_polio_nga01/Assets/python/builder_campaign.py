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

from emod_camp_events import build_node_list, ce_OPV_SIA, ce_OPV_RI, \
                             ce_import_pressure, ce_br_force, \
                             ce_random_numbers, ce_inf_force
from emod_constants import CAMP_FILE

# *****************************************************************************


def campaignBuilder():

    # Variables for this simulation
    START_YEAR = gdata.var_params['start_year']
    RUN_YEARS = gdata.var_params['run_years']

    SIA_COVER = gdata.var_params['sia_base_coverage']
    SIA_RND_SCALE = gdata.var_params['sia_coverage_scale']
    SIA_BASE_TAKE = gdata.var_params['sia_base_vax_take']
    SIA_SCENARIO = gdata.var_params['sia_plan_file']

    RNG_LIST = gdata.var_params['rng_list_offset_yr']
    RNG_VAL = gdata.var_params['rng_list_val']

    RI_START_YR = gdata.var_params['ri_start_yr']

    NODE_DICT = gdata.demog_node

    TIME_MIN = 365.0*(START_YEAR-gdata.base_year)
    TIME_MAX = TIME_MIN + 365.0*RUN_YEARS + gdata.t_step_days

    # Note: campaign module itself is the file object; no Campaign class
    ALL_NODES = gdata.demog_object.node_ids

    # SIA random effects multiplier
    fname = os.path.join('Assets', 'data', 'rand_effect_sia.json')
    with open(fname) as fid01:
        dict_sia_rnd = json.load(fid01)

    for reg_name in dict_sia_rnd:
        p_val = dict_sia_rnd[reg_name]
        nval01 = SIA_RND_SCALE*p_val + np.log(SIA_COVER/(1-SIA_COVER))
        nval02 = 1.0/(1.0+np.exp(-nval01))
        dict_sia_rnd[reg_name] = nval02

    # Apply historic SIA calendar
    fname = os.path.join('Assets', 'data', 'sia_dat.json')
    with open(fname) as fid01:
        dict_sia = json.load(fid01)

    # Apply planned SIA calendar
    if (SIA_SCENARIO):
        fname = os.path.join('Assets', 'data', SIA_SCENARIO+'.json')
        with open(fname) as fid01:
            dict_sia.update(json.load(fid01))

    # Build SIA events
    for sia_name in dict_sia:
        sia_obj = dict_sia[sia_name]

        start_day = sia_obj['date']
        if (start_day < TIME_MIN):
            continue

        if (sia_obj['type'] == 'sabin2'):
            clade = 0
            genome = gdata.boxes_nopv2
            sia_take = SIA_BASE_TAKE
        elif (sia_obj['type'] == 'nopv2'):
            clade = 1
            genome = 0
            sia_take = SIA_BASE_TAKE*gdata.nopv2_sia_take_fac

        for reg_name in dict_sia_rnd:
            cover_val = dict_sia_rnd[reg_name]
            nname_dict = {nname: nname for nname in sia_obj['nodes']}

            n_list_sia = build_node_list([reg_name], nname_dict)
            if (not n_list_sia):
                continue

            n_list = build_node_list(n_list_sia, NODE_DICT)
            camp_event = ce_OPV_SIA(n_list, start_day=start_day,
                                    coverage=cover_val, take=sia_take,
                                    clade=clade, genome=genome)
            camp_module.add(camp_event)

    # Seed infections
    for seed_set in gdata.seed_sets:
        reg_name = seed_set[0]
        seed_time = seed_set[1]
        n_list = build_node_list([reg_name], NODE_DICT)

        # Preserve size of outbreak; select single node for initial location
        n_list = [n_list[-1]]
        cvdpv_gen = gdata.boxes_nopv2+gdata.boxes_sabin2
        start_day = 365.0*(seed_time-gdata.base_year)
        camp_event = ce_import_pressure(n_list, start_day=start_day,
                                        genome=cvdpv_gen,
                                        duration=gdata.seed_inf_dt,
                                        magnitude=gdata.seed_inf_num)
        camp_module.add(camp_event)

    # Time varying birth rate
    start_day = 365.0*(START_YEAR-gdata.base_year)
    for br_tup in gdata.brate_mult_tup:
        BR_MULT_X = br_tup[0]
        BR_MULT_Y = br_tup[1]
        NODE_VALS = br_tup[2]
        camp_event = ce_br_force(NODE_VALS, BR_MULT_X, BR_MULT_Y, start_day)
        camp_module.add(camp_event)

    # Seasonality
    start_day = 365.0*(START_YEAR-gdata.base_year)
    camp_event = ce_inf_force(ALL_NODES, 165.0, 80.0, 1.2, nreps=8,
                              start_day=start_day, dt=gdata.t_step_days)
    camp_module.add(camp_event)

    # Add RI
    ri_age_day = 120.0
    start_day = 365.0*(RI_START_YR-gdata.base_year) - ri_age_day

    with open(os.path.join('Assets', 'data', 'routine.json')) as fid01:
        dict_ri = json.load(fid01)
    if (start_day > TIME_MAX):
        dict_ri = dict()

    for reg_name in dict_ri:
        imm_value = dict_ri[reg_name]
        n_list = build_node_list([reg_name], NODE_DICT)
        camp_event = ce_OPV_RI(n_list, coverage=imm_value, start_day=start_day,
                               base_take=1.0, age_one=ri_age_day, age_std=15.0,
                               clade=1, genome=0)
        camp_module.add(camp_event)

    # Random number stream offset
    for (yr_off, nval) in zip(RNG_LIST, RNG_VAL):
        start_day = 365.0*(START_YEAR-gdata.base_year+yr_off)
        if (nval < 0):
            nval = gdata.sim_index
        camp_event = ce_random_numbers(ALL_NODES, start_day=start_day,
                                       numbers=nval)
        camp_module.add(camp_event)

    # End file construction
    camp_module.save(filename=CAMP_FILE)

    return None

# *****************************************************************************
