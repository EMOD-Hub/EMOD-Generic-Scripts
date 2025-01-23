# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))
from py_assets_common.emod_constants import NUM_SIMS, P_FILE, D_FILE, \
                                            POP_PYR, EXP_C
from py_assets_common.emod_local_proc import pyr_chart

from global_data import base_year

# *****************************************************************************

DIRNAMES = ['experiment_cVDPV2_NGA_100km_baseline']

# *****************************************************************************


def make_fig():

    for dirname in DIRNAMES:

        # Sim outputs
        tpath = os.path.join('..', dirname)

        with open(os.path.join(tpath, D_FILE)) as fid01:
            data_brick = json.load(fid01)

        with open(os.path.join(tpath, P_FILE)) as fid01:
            param_dict = json.load(fid01)

        NSIMS = int(param_dict[NUM_SIMS])
        start_year = param_dict[EXP_C]['start_year']
        run_years = param_dict[EXP_C]['run_years']
        pyr_mat = np.zeros((NSIMS, int(run_years)+1, 20))-1
        year_vec = np.arange(start_year, start_year+run_years+1, dtype=int)
        chart_yrs = year_vec[::6]
        num_charts = chart_yrs.shape[0]

        fig01 = plt.figure(figsize=(8*num_charts, 6))

        data_brick.pop('t_vec')
        data_brick.pop('node_names')
        for sim_idx_str in data_brick:
            sim_idx = int(sim_idx_str)
            pyr_mat[sim_idx, :, :] = np.array(data_brick[sim_idx_str][POP_PYR])

        fidx = (pyr_mat[:, 0, 0] >= 0)

        pyr_mat_avg = np.mean(pyr_mat[fidx, :, :], axis=0)
        pyr_mat_std = np.std(pyr_mat[fidx, :, :], axis=0)

        # Figures - Sims
        for k1 in range(num_charts):

            axs01 = fig01.add_subplot(1, num_charts, k1+1)

            gidx = np.argwhere(year_vec == chart_yrs[k1])[0][0]
            pop_dat = pyr_mat_avg[gidx, :]
            pop_dat_err = pyr_mat_std[gidx, :]

            pyr_chart(axs01, pop_dat, pop_dat_err, year_vec[gidx])

        plt.tight_layout()
        plt.savefig('fig_pyr_{:s}_01.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
