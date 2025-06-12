# *****************************************************************************
#
# *****************************************************************************

import json
import sqlite3

import global_data as gdata

import numpy as np

from emod_postproc_func import post_proc_poppyr, post_proc_cbr, post_proc_R0, \
                               post_proc_sql
from emod_constants import SQL_TIME, SQL_MCW, SQL_AGE, POP_AGE_DAYS, O_FILE, \
                           SQL_FILE

# *****************************************************************************


def application(output_path):

    # Prep output dictionary
    SIM_IDX = gdata.sim_index
    key_str = '{:05d}'.format(SIM_IDX)
    parsed_dat = {key_str: dict()}

    # Sample population pyramid every year
    post_proc_poppyr(output_path, parsed_dat[key_str])

    # Retain annualized count of births
    post_proc_cbr(output_path, parsed_dat[key_str])

    # Estimate R0 from output
    post_proc_R0(output_path, parsed_dat[key_str])

    # Update SQL data
    (dvec_time, _, dvec_mcw, dvec_age) = post_proc_sql(0.0)

    # Yearly timeseries by age
    DAY_BINS = [365]
    START_TIME = 365.0*(gdata.start_year-gdata.base_year)
    BIN_EDGES = np.cumsum(int(gdata.run_years)*DAY_BINS) + START_TIME + 0.5
    BIN_EDGES = np.insert(BIN_EDGES, 0, START_TIME + 0.5)

    inf_dat = np.zeros((int(gdata.run_years), len(POP_AGE_DAYS)-1))
    for k1 in range(inf_dat.shape[1]):
        idx = (dvec_age >= POP_AGE_DAYS[k1]) & \
              (dvec_age < POP_AGE_DAYS[k1+1])
        (inf_yr, tstamps) = np.histogram(dvec_time[idx],
                                         bins=BIN_EDGES,
                                         weights=dvec_mcw[idx])
        inf_dat[:, k1] = inf_yr
    parsed_dat[key_str]['inf_data'] = inf_dat.tolist()

    # Write output dictionary
    with open(O_FILE, 'w') as fid01:
        json.dump(parsed_dat, fid01)

    return None

# *****************************************************************************
