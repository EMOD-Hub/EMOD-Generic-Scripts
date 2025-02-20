# *****************************************************************************
#
# *****************************************************************************

import json
import os
import sqlite3

import global_data as gdata

import numpy as np

from emod_constants import SQL_TIME, SQL_NODE, SQL_AGENT, SQL_LABEL, O_FILE, \
                           SQL_FILE

# *****************************************************************************


def application(output_path):

    # Prep output dictionary
    SIM_IDX = gdata.sim_index
    key_str = '{:05d}'.format(SIM_IDX)
    parsed_dat = {key_str: dict()}

    # Connect to SQL database; retreive new entries
    connection_obj = sqlite3.connect(SQL_FILE)
    cursor_obj = connection_obj.cursor()

    sql_cmd = "SELECT * FROM SIM_EVENTS WHERE SIM_TIME >= {:.1f}".format(0.0)
    cursor_obj.execute(sql_cmd)
    rlist = cursor_obj.fetchall()

    ndata = np.zeros((len(rlist), 4), dtype=int)
    ndata[:, 0] = np.array([val[SQL_TIME] for val in rlist], dtype=int)
    ndata[:, 1] = np.array([val[SQL_NODE] for val in rlist], dtype=int)
    ndata[:, 2] = np.array([val[SQL_AGENT] for val in rlist], dtype=int)
    ndata[:, 3] = np.array([val[SQL_LABEL] for val in rlist], dtype=int)

    # Add tree data to summary output
    parsed_dat[key_str] = ndata.tolist()

    # Write output dictionary
    with open(O_FILE, 'w') as fid01:
        json.dump(parsed_dat, fid01)

    # Validate entries in tree data
    with open(os.path.join(output_path, 'InsetChart.json')) as fid01:
        inset_chart = json.load(fid01)
    ref_inf = np.sum(inset_chart['Channels']['New Infections']['Data'])
    sql_inf = ndata.shape[0]

    # Throw exception if data issues
    if (sql_inf != ref_inf):
        raise Exception('IC = {:d}; SQL = {:d}'.format(ref_inf, sql_inf))

    return None

# *****************************************************************************
