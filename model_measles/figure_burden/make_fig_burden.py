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

from py_assets_common.emod_constants import NUM_SIMS, P_FILE, POP_PYR, D_FILE
from global_data import run_years

# *****************************************************************************

DIRNAMES = ['experiment_sweep_base',
            'experiment_sweep_base_MCV2',
            'experiment_sweep_base_SIAs',]

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
        pyr_mat = np.zeros((NSIMS, int(run_years)+1, 20))-1
        inf_dat = np.zeros((NSIMS, 12*int(run_years)))

        for skey in data_brick:
            if (not skey.isdigit()):
                continue

            sidx = int(skey)
            inf_dat[sidx, :] = np.array(data_brick[skey]['timeseries'])
            pyr_mat[sidx, :, :] = np.array(data_brick[skey][POP_PYR])

        fidx = (pyr_mat[:, 0, 0] >= 0)

        pyr_mat_avg = np.mean(pyr_mat[fidx, :, :], axis=0)
        tpop_avg = np.sum(pyr_mat_avg, axis=1)
        tpop_xval = np.arange(len(tpop_avg))

        # Figures
        fig01 = plt.figure(figsize=(8, 6))

        # Figures - Sims - Infections
        axs01 = fig01.add_subplot(1, 1, 1)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)

        xval = np.arange(0, run_years, 1/12) + 1/24
        pops = np.interp(xval, tpop_xval, tpop_avg)
        inf_dat = inf_dat[fidx, :]/pops*1e5
        yval = np.mean(inf_dat, axis=0)

        infDatSetSort = np.sort(inf_dat, axis=0)
        infDatSetSort = infDatSetSort

        for patwid in [0.45, 0.375, 0.25]:
            xydat = np.zeros((2*infDatSetSort.shape[1], 2))
            xydat[:, 0] = np.hstack((xval, xval[::-1]))
            tidx = int((0.5-patwid)*infDatSetSort.shape[0])
            xydat[:, 1] = np.hstack((infDatSetSort[tidx, :],
                                     infDatSetSort[-tidx, ::-1]))

            polyShp = patch.Polygon(xydat, facecolor='C0',
                                    alpha=0.7-patwid, edgecolor=None)
            axs01.add_patch(polyShp)

        axs01.plot(xval, yval, color='C0', linewidth=2)
        axs01.set_ylabel('Monthly Infections per 100k', fontsize=16)
        axs01.set_xlabel('Year', fontsize=16)

        plt.tight_layout()
        plt.savefig('fig_clouds_{:s}01.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
