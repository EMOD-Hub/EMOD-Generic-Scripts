# *****************************************************************************
#
# *****************************************************************************

import json

import global_data as gdata

import numpy as np

import emod_api.campaign as camp_module

from emod_camp_events import ce_SIA
from emod_postproc_func import post_proc_sql

# *****************************************************************************


def application(timestep):

    TIME_VAL = float(timestep)

    PREV_TIME = gdata.prev_proc_time
    NODE_IDS = gdata.node_idval
    CASE_THRESH = gdata.var_params['adm01_case_threshold']
    LOG_REP_RATE = gdata.var_params['log10_min_reporting']

    ADM01_DICT = gdata.adm01_idlist
    ADM01_CASE = gdata.adm01_cases

    # First call to in-process; initialize as needed
    if (gdata.first_call_bool):
        gdata.first_call_bool = False

        gdata.data_vec_time = np.array(list())
        gdata.data_vec_node = np.array(list())
        gdata.data_vec_mcw = np.array(list())

        gdata.adm01_cases = {adm01: 0 for adm01 in ADM01_DICT}

        adm02_rate = {adm02: np.random.beta(0.33, 8.0) 
                      for adm02 in gdata.adm02_idlist}

        with open('adm02_obsrate.json', 'w') as fid01:
            json.dump(adm02_rate, fid01)

        max_node_id = np.max(gdata.demog_object.node_ids)
        nobs_list = np.zeros(max_node_id, dtype=float)
        for nname in NODE_IDS:
            for adm02_name in adm02_rate:
                if (nname.startswith(adm02_name+':')):
                    nobs_list[NODE_IDS[nname]-1] = adm02_rate[adm02_name]

        gdata.nobs_vec = nobs_list

        print("Hello from in-process at time {:.1f}".format(TIME_VAL))

    # Only evaluate every month-ish
    if (((TIME_VAL % 365.0) % 30.0) < 29.0):
        return None

    # Update SQL data
    (dvec_time, dvec_node, dvec_mcw, _) = post_proc_sql(PREV_TIME)

    gdata.data_vec_time = np.append(gdata.data_vec_time, dvec_time)
    gdata.data_vec_node = np.append(gdata.data_vec_node, dvec_node)
    gdata.data_vec_mcw = np.append(gdata.data_vec_mcw, dvec_mcw)

    # Record timestamp of previous process
    gdata.prev_proc_time = TIME_VAL

    # Only do response post-2020
    obr_time = 365.0*(gdata.start_year_obr - gdata.base_year)
    if (TIME_VAL > obr_time):
        for k1 in range(gdata.nobs_vec.shape[0]):
            nid_val = k1 + 1
            gidx = (dvec_node == nid_val)
            min_rate = np.power(10.0, LOG_REP_RATE)
            nreprate = np.maximum(gdata.nobs_vec[nid_val-1], min_rate)
            obs_inf = np.sum(dvec_mcw[gidx])*nreprate
            for adm01 in ADM01_DICT:
                if (nid_val in ADM01_DICT[adm01]):
                    gdata.adm01_cases += obs_inf

    # Reactive campaign
    targ_nodes = list()
    for adm01 in ADM01_DICT:
        if (gdata.adm01_cases[adm01] > CASE_THRESH):
            targ_nodes.extend(ADM01_DICT[adm01])
            gdata.adm01_cases = 0

    # New campaign file
    camp_module.reset()
    CAMP_FILE = 'campaign_{:05d}.json'.format(int(TIME_VAL))
    ALL_NODES = gdata.demog_object.node_ids

    if (targ_nodes):
        sia_day = TIME_VAL+30.0
        camp_event = ce_SIA(ALL_NODES, start_day=sia_day, coverage=0.50,
                            yrs_min=0.75, yrs_max=5.0)
        camp_module.add(camp_event)

        camp_module.save(filename=CAMP_FILE)
        return CAMP_FILE

    return None

# *****************************************************************************
