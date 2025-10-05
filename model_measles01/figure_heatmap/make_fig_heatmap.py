# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))

from py_assets_common.emod_constants import EXP_V, NUM_SIMS, P_FILE, POP_PYR, \
                                            D_FILE
from global_data import run_years, AGE_HIST_BINS, IHME_MORT_X, IHME_MORT_Y
                        

# *****************************************************************************

DIRNAMES = ['experiment_surf_base',
            'experiment_surf_ideal']

# *****************************************************************************


def make_fig():

    for dirname in DIRNAMES:

        # Sim outputs
        tpath = os.path.join('..', dirname)

        with open(os.path.join(tpath, D_FILE)) as fid01:
            data_brick = json.load(fid01)

        with open(os.path.join(tpath, P_FILE)) as fid01:
            param_dict = json.load(fid01)

        nsims = int(param_dict[NUM_SIMS])

        inf_dat = np.zeros((nsims, 12*int(run_years)))
        age_dat = np.zeros((nsims, int(run_years), len(AGE_HIST_BINS)-1))
        mort_dat = np.zeros((nsims, int(run_years)))
        pyr_mat = np.zeros((nsims, int(run_years)+1, 20))-1

        mcv1_vec = np.array(param_dict[EXP_V]['MCV1'])
        mcv1_lvl = np.unique(mcv1_vec)

        mcv1_age_vec = np.array(param_dict[EXP_V]['MCV1_age'])
        mcv1_age_lvl = np.unique(mcv1_age_vec).tolist()

        xyrs = np.arange(0, run_years, 1) + 1/2
        xages = AGE_HIST_BINS[1:] + np.diff(AGE_HIST_BINS)/2
        mort_prob = np.interp(xages, IHME_MORT_X, IHME_MORT_Y)

        for skey in data_brick:
            if (not skey.isdigit()):
                continue

            sidx = int(skey)
            inf_dat[sidx, :] = np.array(data_brick[skey]['timeseries'])
            age_dat[sidx, :, :] = np.array(data_brick[skey]['age_data'])
            pyr_mat[sidx, :, :] = np.array(data_brick[skey][POP_PYR])
            mort_dat[sidx, :] = np.sum(age_dat[sidx, :, :]*mort_prob, axis=1)

        fidx = (pyr_mat[:, 0, 0] >= 0)

        pyr_mat_avg = np.mean(pyr_mat[fidx, :, :], axis=0)
        tpop_avg = np.sum(pyr_mat_avg, axis=1)
        tpop_xval = np.arange(len(tpop_avg))

        pops = np.interp(xyrs, tpop_xval, tpop_avg)
        mort_nrm = mort_dat[fidx, :]/pops*1e5

        # Subset to only valid simulations
        mcv1_vec = mcv1_vec[fidx]
        mcv1_age_vec = mcv1_age_vec[fidx]

        # Figures
        fig01 = plt.figure(figsize=(10, 6))

        # Figures - Sims - Infections
        axs01 = fig01.add_subplot(1, 1, 1)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)

        # Binning
        cval = np.mean(mort_nrm[:, -10:], axis=1)/12.0
        (xvec, yvec) = np.meshgrid(mcv1_lvl, mcv1_age_lvl)
        zmat = np.zeros(xvec.shape)
        for k1 in range(xvec.shape[0]):
            for k2 in range(xvec.shape[1]):
                gidx = (mcv1_vec == xvec[k1, k2]) & \
                       (mcv1_age_vec == yvec[k1, k2])
                zmat[k1, k2] = np.mean(cval[gidx])

        plt.contour(12*yvec/365.0, xvec, zmat,
                    levels=[0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0],
                    linewidths=3, vmin=0, vmax=3)
        axs01.scatter(12*mcv1_age_vec/365.0, mcv1_vec, c=cval, vmin=0, vmax=3)

        virmap = plt.get_cmap('viridis_r')
        cbar_handle = plt.colorbar(cm.ScalarMappable(cmap=virmap),
                                   ax=axs01, shrink=0.75)

        ticloc = [0, 1/6, 2/6, 3/6, 4/6, 5/6, 1]
        ticlab = ['3.0', '2.5', '2.0', '1.5', '1.0', '0.5', '0.0']

        cbar_handle.set_ticks(ticks=ticloc)
        cbar_handle.set_ticklabels(ticlab)
        cbar_handle.set_label('Monthly Mortality per-100k',
                              fontsize=14, labelpad=10)
        cbar_handle.ax.tick_params(labelsize=14)

        axs01.set_xlabel('MCV1 Age Policy (months)', fontsize=16)
        axs01.set_ylabel('MCV Coverage', fontsize=16)
        axs01.set_xlim(4, 18)
        axs01.set_ylim(0.2, 1.0)
        axs01.tick_params(axis='both', which='major', labelsize=14)

        plt.tight_layout()
        plt.savefig('fig_heatmap_{:s}01.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
