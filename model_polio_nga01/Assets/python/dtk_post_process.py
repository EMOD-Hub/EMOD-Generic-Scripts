# *****************************************************************************
#
# *****************************************************************************

import json
import os

import global_data as gdata

import numpy as np

from emod_postproc_func import post_proc_poppyr
from emod_constants import RST_FILE, RST_TIME, RST_NODE, RST_CLADE, \
                           RST_GENOME, RST_NEW_INF

# *****************************************************************************


def application(output_path):

    # Variables for this simulation
    START_YEAR = gdata.var_params['start_year']
    RUN_YEARS = gdata.var_params['run_years']
    cVDPV_genome = gdata.boxes_sabin2 + gdata.boxes_nopv2

    ADM02_DICT = gdata.adm02_idlist
    ADM02_LIST = sorted(list(ADM02_DICT.keys()))
    ADM02_IDX = {val[0]: val[1] for val in
                 zip(ADM02_LIST, range(len(ADM02_LIST)))}

    # Prep output dictionary
    SIM_IDX = gdata.sim_index
    key_str = '{:05d}'.format(SIM_IDX)
    parsed_dat = {key_str: dict()}

    # Sample population pyramid every year
    post_proc_poppyr(output_path, parsed_dat[key_str])

    # Timesteps
    t_ini = 365.0*(START_YEAR-gdata.base_year)
    t_end = t_ini + 365.0*RUN_YEARS
    t_vec = np.arange(t_ini, t_end+0.5, gdata.t_step_days)

    # Post-process strain reporter
    infdat = np.loadtxt(os.path.join(output_path, RST_FILE),
                        delimiter=',', skiprows=1, ndmin=2)

    # Construct csv file for cVDPV infections
    dbrick0 = np.zeros((len(ADM02_LIST), t_vec.shape[0]))

    if (infdat.shape[0] > 0):
        for adm02 in ADM02_LIST:
            odex = ADM02_IDX[adm02]
            gidx = (infdat[:, RST_CLADE] == 0)
            gidx = (infdat[:, RST_GENOME] == cVDPV_genome) & gidx
            gidx = np.isin(infdat[:, RST_NODE], ADM02_DICT[adm02]) & gidx
            subdat = infdat[gidx, :]
            for k1 in range(subdat.shape[0]):
                tdex = int((subdat[k1, RST_TIME]-t_ini)/gdata.t_step_days)
                dbrick0[odex, tdex] += subdat[k1, RST_NEW_INF]

    # Log data for local machine
    parsed_dat[key_str]['infmat'] = dbrick0.tolist()
    parsed_dat['t_vec'] = t_vec.tolist()
    parsed_dat['node_names'] = ADM02_IDX

    # Blank output for sims without outbreak
    tot_inf = np.sum(dbrick0, axis=0)
    cum_inf = np.cumsum(tot_inf)
    if (cum_inf[-1] < gdata.init_ob_thresh):
        parsed_dat = dict()

    if (RUN_YEARS >= 3 and cum_inf[-1] < 80e3):
        parsed_dat = dict()

    if (RUN_YEARS >= 5 and cum_inf[-1] < 500e3):
        parsed_dat = dict()

    # Write output dictionary
    with open('parsed_out.json', 'w') as fid01:
        json.dump(parsed_dat, fid01)

    return None

# *****************************************************************************
