# *****************************************************************************
#
# *****************************************************************************

import json

import global_data as gdata

import numpy as np

import emod_api.campaign as camp_module

from emod_camp_events import ce_SIA
from emod_post_proc_func import post_proc_sql

# *****************************************************************************


def application(timestep):

    TIME_VAL = float(timestep)

    PREV_TIME = gdata.prev_proc_time
    NODE_DICT = gdata.demog_node
    CASE_THRESH = gdata.var_params['adm01_case_threshold']
    LOG_REP_RATE = gdata.var_params['log10_min_reporting']

    # First call to in-process; initialize as needed
    if (gdata.first_call_bool):
        gdata.first_call_bool = False
        gdata.data_vec_time = np.array(list())
        gdata.data_vec_node = np.array(list())
        gdata.data_vec_mcw = np.array(list())

        adm01_name = list(set([':'.join(val.split(':')[:3]) for val in NODE_DICT]))
        adm01_list = [[val, list(), 0] for val in adm01_name]
        for nname in NODE_DICT:
            for adm01_tup in adm01_list:
                if (nname.startswith(adm01_tup[0]+':')):
                    adm01_tup[1].append(NODE_DICT[nname])

        adm02_name = list(set([':'.join(val.split(':')[:4]) for val in NODE_DICT]))
        adm02_rate = {val: np.random.beta(0.33, 8.0) for val in adm02_name}

        with open('adm02_obsrate.json', 'w') as fid01:
            json.dump(adm02_rate, fid01)

        nobs_list = np.zeros(gdata.max_node_id, dtype=float)
        for nname in NODE_DICT:
            for adm02_name in adm02_rate:
                if (nname.startswith(adm02_name+':')):
                    nobs_list[NODE_DICT[nname]-1] = adm02_rate[adm02_name]

        gdata.adm01_list = adm01_list
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

    # Only do response starting in 2020
    if (TIME_VAL > 43800.0):
        for k1 in range(gdata.max_node_id):
            nid_val = k1 + 1
            gidx = (dvec_node == nid_val)
            min_rate = np.power(10.0, LOG_REP_RATE)
            nreprate = np.maximum(gdata.nobs_vec[nid_val-1], min_rate)
            obs_inf = np.sum(dvec_mcw[gidx])*nreprate
            for adm01_tup in gdata.adm01_list:
                if (nid_val in adm01_tup[1]):
                    adm01_tup[2] += obs_inf

    # Reactive campaign
    targ_nodes = list()
    for adm01_tup in gdata.adm01_list:
        if (adm01_tup[2] > CASE_THRESH):
            targ_nodes.extend(adm01_tup[1])
            adm01_tup[2] = 0

    # New campaign file
    camp_module.reset()
    CAMP_FILE = 'campaign_{:05d}.json'.format(int(TIME_VAL))
    ALL_NODES = gdata.demog_object.node_ids

    if (targ_nodes):
        sia_day = TIME_VAL+30.0
        camp_event = ce_SIA(ALL_NODES, start_day=sia_day, coverage=0.50,
                            yrs_min=0.75, yrs_max=5.0)
        camp_module.add(camp_event)

    if (camp_module.campaign_dict["Events"]):
        camp_module.save(filename=CAMP_FILE)
        return CAMP_FILE

    return None

# *****************************************************************************
