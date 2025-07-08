# *****************************************************************************
#
# *****************************************************************************

import os
import json

import global_data as gdata

import numpy as np

from emod_analysis import norpois_opt
from emod_postproc_func import post_proc_poppyr, post_proc_sql
from emod_constants import O_FILE, MO_DAYS

# *****************************************************************************


def application(output_path):

    # Variables for this simulation
    END_YEAR = gdata.var_params['end_year']

    # Prep output dictionary
    SIM_IDX = gdata.sim_index
    key_str = '{:05d}'.format(SIM_IDX)
    parsed_dat = {key_str: dict()}

    # Sample population pyramid every year
    post_proc_poppyr(output_path, parsed_dat[key_str])

    # Update SQL data
    PREV_TIME = gdata.prev_proc_time
    (dvec_time, dvec_node, dvec_mcw, _) = post_proc_sql(PREV_TIME)

    gdata.data_vec_time = np.append(gdata.data_vec_time, dvec_time)
    gdata.data_vec_node = np.append(gdata.data_vec_node, dvec_node)
    gdata.data_vec_mcw = np.append(gdata.data_vec_mcw, dvec_mcw)

    # Aggregate new infections by month
    START_TIME = 365.0*(gdata.start_year_log-gdata.base_year)
    NUM_YEARS = int(END_YEAR - gdata.start_year_log)
    BIN_EDGES = np.cumsum(NUM_YEARS*MO_DAYS) + START_TIME + 0.5
    BIN_EDGES = np.insert(BIN_EDGES, 0, START_TIME + 0.5)

    (inf_mo, tstamps) = np.histogram(gdata.data_vec_time,
                                     bins = BIN_EDGES,
                                     weights = gdata.data_vec_mcw)

    # Monthly timeseries
    parsed_dat[key_str]['timeseries'] = inf_mo.tolist()

    # Common output data
    parsed_dat['tstamps'] = (np.diff(tstamps)/2.0 + tstamps[:-1]).tolist()

    # Reference data
    tpath = os.path.join('Assets', 'data', 'GHA_epi.json')
    with open(tpath) as fid01:
        ref_dat = json.load(fid01)
    ref_cases = np.array(ref_dat['cases_mo'])
    ref_start = ref_dat['start_year']
    ref_years = int(np.ceil(ref_cases.shape[0]/12))

    # Aggregate new infections by month
    START_TIME = 365.0*(ref_start-gdata.base_year)
    BIN_EDGES = np.cumsum(ref_years*MO_DAYS) + START_TIME + 0.5
    BIN_EDGES = np.insert(BIN_EDGES, 0, START_TIME + 0.5)

    (inf_mo, _) = np.histogram(gdata.data_vec_time,
                               bins = BIN_EDGES,
                               weights = gdata.data_vec_mcw)

    # Calibration score from timeseries data
    (obj_val, scal_vec) = norpois_opt(ref_cases, inf_mo)
    parsed_dat[key_str]['cal_val'] = float(obj_val)
    parsed_dat[key_str]['rep_rate'] = float(scal_vec)

    # Write output dictionary
    with open(O_FILE, 'w') as fid01:
        json.dump(parsed_dat, fid01)

    return None

# *****************************************************************************
