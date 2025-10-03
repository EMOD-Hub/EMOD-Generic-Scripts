# *****************************************************************************

import json
import os
import sys

import numpy as np

from aux_defaultgrid import svg_defaultgrid

# Ought to go in emodpy
LOCAL_PATH = os.path.abspath(os.path.join('..', '..', 'local_python'))
sys.path.insert(0, LOCAL_PATH)
from py_assets_common.emod_constants import EXP_V, P_FILE

# *****************************************************************************

DIRNAME = 'experiment_network01'

# *****************************************************************************


def make_fig():

    # Sim outputs
    tpath = os.path.join('..', DIRNAME)

    with open(os.path.join(tpath, 'data_brick.json')) as fid01:
        data_brick = json.load(fid01)

    with open(os.path.join(tpath, P_FILE)) as fid01:
        param_dict = json.load(fid01)

    net_coef = np.array(param_dict[EXP_V]['net_coef'])
    netc_vals = np.unique(net_coef)

    # One animation for each coefficient value
    for k1 in range(netc_vals.shape[0]):
        sim_idx = np.argwhere(net_coef == netc_vals[k1])[0][0]
        inf_dat = data_brick['{:05d}'.format(sim_idx)]['inf_data']
        f_lab = '{:1.0e}'.format(netc_vals[k1])
        f_lab = f_lab.replace('+', '')
        svg_defaultgrid(inf_dat, f_lab)

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
