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
                             ce_random_numbers, ce_inf_force, ce_inf_mod
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
    SIA_END_YR = gdata.var_params['sia_end_yr']

    RNG_LIST = gdata.var_params['rng_list_offset_yr']
    RNG_VAL = gdata.var_params['rng_list_val']

    RI_START_YR = gdata.var_params['ri_start_yr']

    R0_MOD = gdata.var_params['R0_mods']

    NODE_DICT = gdata.node_idval

    DAY_MIN = 365.0*(START_YEAR-gdata.base_year)
    DAY_MAX = DAY_MIN + 365.0*RUN_YEARS + gdata.t_step_days
    DAY_MAX_SIA = 365.0*(SIA_END_YR-gdata.base_year)

    # Note: campaign module itself is the file object; no Campaign class
    ALL_NODES = gdata.demog_object.node_ids

    # SIA random effects value
    fname = os.path.join('Assets', 'data', 'sia_rnd_effect.json')
    with open(fname) as fid01:
        dict_sia_rnd = json.load(fid01)

    # SIA coverage-by-node {node_id: sia_coverage}
    sia_cover_dict = {NODE_DICT[nname]: SIA_COVER for nname in NODE_DICT}

    # Update SIA coverage-by-node with random effects value
    for reg_name in dict_sia_rnd:
        eff_val = dict_sia_rnd[reg_name]
        int_val = SIA_RND_SCALE*eff_val + np.log(SIA_COVER/(1-SIA_COVER))
        cvr_val = 1.0/(1.0+np.exp(-int_val))

        n_list = build_node_list([reg_name], NODE_DICT)
        sia_cover_dict.update({n_id: cvr_val for n_id in n_list})
    gdata.sia_cover_dict = sia_cover_dict

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
        sia_day = dict_sia[sia_name]['date']
        if (sia_day < DAY_MIN or sia_day > DAY_MAX or sia_day > DAY_MAX_SIA):
            continue

        if (dict_sia[sia_name]['type'] == 'sabin2'):
            clade = 0
            genome = gdata.boxes_nopv2
            sia_take = SIA_BASE_TAKE
        elif (dict_sia[sia_name]['type'] == 'nopv2'):
            clade = 1
            genome = 0
            sia_take = SIA_BASE_TAKE*gdata.nopv2_sia_take_fac

        age_yr_max = dict_sia[sia_name]['age_yr_max']

        n_list = build_node_list(dict_sia[sia_name]['nodes'], NODE_DICT)
        if (n_list):
            n_dict = {nid: sia_cover_dict[nid] for nid in n_list}
            camp_event = ce_OPV_SIA(n_dict, start_day=sia_day, take=sia_take,
                                    yrs_min=0.2, yrs_max=age_yr_max,
                                    clade=clade, genome=genome)
            camp_module.add(camp_event)

    # Seed infections
    for seed_set in gdata.seed_sets:
        cvdpv_gen = gdata.boxes_nopv2+gdata.boxes_sabin2
        reg_name = seed_set[0]
        seed_time = seed_set[1]

        n_list = build_node_list([reg_name], NODE_DICT)
        if (n_list):
            n_list = [n_list[-1]]  # Select single node for initial location
            start_day = 365.0*(seed_time-gdata.base_year)
            camp_event = ce_import_pressure(n_list, start_day=start_day,
                                            genome=cvdpv_gen,
                                            duration=gdata.seed_inf_dt,
                                            magnitude=gdata.seed_inf_num)
            camp_module.add(camp_event)

    # Time varying birth rate
    start_day = 365.0*(START_YEAR-gdata.base_year)
    for brate_mult_tup in gdata.brate_mult_tup_list:
        BR_MULT_X = brate_mult_tup[0]
        BR_MULT_Y = brate_mult_tup[1]
        NODE_VALS = brate_mult_tup[2]
        camp_event = ce_br_force(NODE_VALS, BR_MULT_X, BR_MULT_Y, start_day)
        camp_module.add(camp_event)

    # Add seasonality
    start_day = 365.0*(START_YEAR-gdata.base_year)

    fname = os.path.join('Assets', 'data', 'r0_seasonality.json')
    with open(fname) as fid01:
        dict_r0_season = json.load(fid01)

    for reg_name in dict_r0_season:
        n_list = build_node_list([reg_name], NODE_DICT)
        if (n_list):
            camp_event = ce_inf_force(n_list,
                                      dict_r0_season[reg_name]['day_start'],
                                      dict_r0_season[reg_name]['day_length'],
                                      dict_r0_season[reg_name]['mult_val'],
                                      start_day=start_day,
                                      dt=gdata.t_step_days)
            camp_module.add(camp_event)

    # Ad hoc R0 modifiers
    for mod_tuple in R0_MOD:
        n_list = build_node_list(mod_tuple[0], NODE_DICT)
        start_day = 365.0*(mod_tuple[1]-gdata.base_year)
        dt_days = 365.0*mod_tuple[2]
        mult_val = mod_tuple[3]
        camp_event = ce_inf_mod(n_list, start_day=start_day,
                                dt_days=dt_days, mult_val=mult_val)
        camp_module.add(camp_event)

    # Add RI
    ri_age_day = 90.0
    start_day = 365.0*(RI_START_YR-gdata.base_year) - ri_age_day

    with open(os.path.join('Assets', 'data', 'routine_dat.json')) as fid01:
        dict_ri = json.load(fid01)
    if (start_day > DAY_MAX):
        dict_ri = dict()

    for reg_name in dict_ri:
        imm_value = dict_ri[reg_name]
        n_list = build_node_list([reg_name], NODE_DICT)
        if (n_list):
            camp_event = ce_OPV_RI(n_list, coverage=imm_value,
                                   start_day=start_day, base_take=1.0,
                                   age_one=ri_age_day, age_std=15.0,
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
