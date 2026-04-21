# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))

from py_assets_common.emod_constants import NUM_SIMS, P_FILE, POP_PYR, \
                                            D_FILE, EXP_C
from global_data import run_years, t_step_days

# *****************************************************************************

DIRNAMES = ['experiment_none_ri70_nodes002',
            'experiment_none_ri70_nodes020']

# *****************************************************************************


def make_fig():

    # Figures
    fig01 = plt.figure(figsize=(8, 6))

    # Figures - Sims - Prevalance
    axs01 = fig01.add_subplot(1, 1, 1)
    plt.sca(axs01)

    axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
    axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
    axs01.set_axisbelow(True)

    for dirname in DIRNAMES:

        # Sim outputs
        tpath = os.path.join('..', dirname)

        with open(os.path.join(tpath, D_FILE)) as fid01:
            data_brick = json.load(fid01)

        with open(os.path.join(tpath, P_FILE)) as fid01:
            param_dict = json.load(fid01)

        NSIMS = int(param_dict[NUM_SIMS])
        NNODES = int(param_dict[EXP_C]["num_nodes"])
        tvec = np.arange(0, 365*run_years, t_step_days)

        inf_dat = np.zeros((NSIMS, NNODES, len(tvec)))
        pyr_mat = np.zeros((NSIMS, int(run_years)+1, 20))-1

        for skey in data_brick:
            if (not skey.isdigit()):
                continue

            sim_idx = int(skey)

            inf_col = np.array(data_brick[skey]['rst'])  # (Time, Node, Infections)
            for k1 in range(inf_col.shape[0]):
                time_idx = int(inf_col[k1, 0]/t_step_days)
                node_idx = int(inf_col[k1, 1] - 1)  # Node indicies are 1-based
                inf_dat[sim_idx, node_idx, time_idx] = inf_col[k1, 2]
            pyr_mat[sim_idx, :, :] = np.array(data_brick[skey][POP_PYR])

        fidx = (pyr_mat[:, 0, 0] >= 0)

        non_zero_dat = (inf_dat[fidx, :, :] > 0)
        non_zero_avg = np.mean(non_zero_dat, axis=(0,1))  # Lump sims and nodes
        non_zero_std = np.std(non_zero_dat, axis=(0,1))  # Lump sims and nodes

        axs01.plot(tvec/365.0, non_zero_avg, linewidth=2,
                   label='Nodes {:d}'.format(NNODES))
        axs01.set_ylabel('Non-zero prevalence (fraction)', fontsize=16)
        axs01.set_xlabel('Year', fontsize=16)

    plt.tight_layout()
    plt.legend()
    plt.savefig('fig_prevalence_01.png'.format(dirname))
    plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
