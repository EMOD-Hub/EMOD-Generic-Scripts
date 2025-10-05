# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))

from py_assets_common.emod_constants import EXP_C, NUM_SIMS, P_FILE, POP_PYR, \
                                            D_FILE
from py_assets_common.emod_local_proc import pyr_chart

# *****************************************************************************

DIRNAMES = ['experiment_estimates01',
            'experiment_projections01']

# *****************************************************************************


def make_fig():

    # Sim outputs
    for dirname in DIRNAMES:
        tpath = os.path.join('..', dirname)

        with open(os.path.join(tpath, D_FILE)) as fid01:
            data_brick = json.load(fid01)

        with open(os.path.join(tpath, P_FILE)) as fid01:
            param_dict = json.load(fid01)

        nsims = int(param_dict[NUM_SIMS])
        init_year = int(param_dict[EXP_C]['start_year'])
        num_years = int(param_dict[EXP_C]['num_years'])
        pop_dat_str = param_dict[EXP_C]['pop_dat_file']

        pyr_mat = np.zeros((nsims, num_years+1, 20))-1
        year_vec = np.arange(init_year, init_year+num_years+1)
        chart_yrs = year_vec[np.mod(year_vec, 10) == 0]
        num_charts = chart_yrs.shape[0]

        fig01 = plt.figure(figsize=(6*num_charts, 10))

        for sim_idx_str in data_brick:
            sim_idx = int(sim_idx_str)
            pyr_mat[sim_idx, :, :] = np.array(data_brick[sim_idx_str][POP_PYR])

        fidx = (pyr_mat[:, 0, 0] >= 0)

        popfile = 'pop_dat_{:s}.csv'.format(pop_dat_str)
        fname_pop = os.path.join('..', 'Assets', 'data', popfile)
        pop_input = np.loadtxt(fname_pop, dtype=int, delimiter=',')

        year_vec_dat = pop_input[0, :]
        pop_mat_dat = pop_input[1:, :]

        pyr_mat_avg = np.mean(pyr_mat[fidx, :, :], axis=0)
        pyr_mat_std = np.std(pyr_mat[fidx, :, :], axis=0)

        # Figures - Sims
        for k1 in range(num_charts):

            gidx = np.argwhere(year_vec == chart_yrs[k1])[0][0]
            pop_dat = pyr_mat_avg[gidx, :]
            pop_dat_err = pyr_mat_std[gidx, :]

            axs01 = fig01.add_subplot(2, num_charts, k1+1)
            plt.sca(axs01)
            pyr_chart(axs01, pop_dat, pop_dat_err, chart_yrs[k1])

            if (k1 == num_charts-1):
                axs02 = axs01.twinx()
                axs02.set_ylabel('Simulation', fontsize=24)
                axs02.set_yticks(ticks=[0, 1])
                axs02.set_yticklabels(['', ''])

        # Figures - Reference
        for k1 in range(num_charts):

            gidx = np.argwhere(year_vec_dat == chart_yrs[k1])[0][0]
            pop_dat = pop_mat_dat[:-1, gidx]
            pop_dat_err = 0*pop_dat

            axs01 = fig01.add_subplot(2, num_charts, k1+1+num_charts)
            plt.sca(axs01)
            pyr_chart(axs01, pop_dat, pop_dat_err, chart_yrs[k1])

            if (k1 == num_charts-1):
                axs02 = axs01.twinx()
                axs02.set_ylabel('Reference', fontsize=24)
                axs02.set_yticks(ticks=[0, 1])
                axs02.set_yticklabels(['', ''])

        # Save figure
        plt.tight_layout()
        plt.savefig('fig_pyr_{:s}01_{:s}.png'.format(pop_dat_str, dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
