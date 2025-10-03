# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'local_python')))
from py_assets_common.emod_constants import EXP_V, NUM_SIMS, P_FILE, D_FILE

# *****************************************************************************

DIRNAME = 'experiment_transtree01'

# *****************************************************************************


def rec_tree(tree_dat, targ_id):

    subdat = tree_dat[tree_dat[:, 3] == targ_id, :]
    subtree = list()

    for k1 in range(subdat.shape[0]):
        addleaf = rec_tree(tree_dat, subdat[k1, 2])
        subtree.append([subdat[k1, 0], addleaf])

    return subtree

# *****************************************************************************


def print_branch(axs_id, leaf_list, c_yval, pxy=None):

    sort_leaf = sorted(leaf_list, key=lambda val: val[0], reverse=True)

    for leaf_val in sort_leaf:
        axs_id.plot(leaf_val[0], c_yval, marker='.', lw=0, c='r', ms=0.5)
        new_pxy = (leaf_val[0], c_yval)
        if (pxy):
            if (pxy[0] < new_pxy[0]):
                xdat = [pxy[0]+0.25, pxy[0]+0.5, pxy[0]+0.5, new_pxy[0]-0.25]
                ydat = [pxy[1],      pxy[1],     new_pxy[1], new_pxy[1]]
            else:  # Secondary infection on same timestep; edge case in seeding
                xdat = [pxy[0],   new_pxy[0]]
                ydat = [pxy[1]+1, new_pxy[1]-1]
            axs_id.plot(xdat, ydat, lw=0.5, c='k')
        c_yval = print_branch(axs_id, leaf_val[1], c_yval+5, new_pxy)

    return c_yval

# *****************************************************************************


def make_fig():

    # Sim outputs
    tpath = os.path.join('..', DIRNAME)

    with open(os.path.join(tpath, D_FILE)) as fid01:
        data_brick = json.load(fid01)

    with open(os.path.join(tpath, P_FILE)) as fid01:
        param_dict = json.load(fid01)

    BI_VAR = np.array(param_dict[EXP_V]['base_inf_stddev_mult'])
    nsims = int(param_dict[NUM_SIMS])

    tree_list = [np.array([]) for k1 in range(nsims)]
    for idx_str in data_brick:
        idx_val = int(idx_str)
        tree_list[idx_val] = np.array(data_brick[idx_str])

    # Simulations with outbreak
    obr_bool = np.array([tree.shape[0] > 100 for tree in tree_list])

    # Number of simulations at each level
    sim_num = np.array([np.sum(BI_VAR == 0.0),
                        np.sum(BI_VAR == 1.0)])

    # Number of simulations with outbreak at each level
    obr_num = np.array([np.sum(obr_bool[BI_VAR == 0.0]),
                        np.sum(obr_bool[BI_VAR == 1.0])])

    print(np.around(obr_num/sim_num, 3))

    gdex0 = np.argwhere(obr_bool & (BI_VAR == 0.0))[0, 0]
    gdex1 = np.argwhere(obr_bool & (BI_VAR == 1.0))[0, 0]

    # Figure setup
    fig01 = plt.figure(figsize=(48, 36))

    axs01 = fig01.add_subplot(1, 2, 1, label=None)
    plt.sca(axs01)
    axs01.tick_params(labelcolor='w',
                      top=False, bottom=False,
                      left=False, right=False)

    tree_dat = tree_list[gdex0]
    x_max = np.max(tree_dat[:, 0])
    axs01.set_xlabel('Time', fontsize=14)
    axs01.set_title('Constant Rate\nSecondary Infection', fontsize=16)
    axs01.set_xlim(-10, x_max + 10)

    print_branch(axs01, rec_tree(tree_dat, 0), 0)

    axs01 = fig01.add_subplot(1, 2, 2, label=None)
    plt.sca(axs01)
    axs01.tick_params(labelcolor='w',
                      top=False, bottom=False,
                      left=False, right=False)

    tree_dat = tree_list[gdex1]
    x_max = np.max(tree_dat[:, 0])
    axs01.set_xlabel('Time', fontsize=14)
    axs01.set_title('Exponential Rate\nSecondary Infection', fontsize=16)
    axs01.set_xlim(-10, x_max + 10)

    print_branch(axs01, rec_tree(tree_dat, 0), 0)

    # Generate figures
    plt.tight_layout()
    plt.savefig('fig_transtree01.png')
    plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
