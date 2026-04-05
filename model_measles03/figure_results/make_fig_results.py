# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))

from py_assets_common.emod_constants import EXP_C, EXP_V, NUM_SIMS, \
                                            P_FILE, POP_PYR, D_FILE
from global_data import run_years, AGE_HIST_BINS, IHME_MORT_X, IHME_MORT_Y

# *****************************************************************************

DIRNAMES = ['experiment_none_ri70',
            'experiment_sia3yr60cv_ri70']

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
        dose_dat = np.zeros((nsims, int(run_years)))
        pyr_mat = np.zeros((nsims, int(run_years)+1, 20))-1

        mcv1_val = np.array(param_dict[EXP_C]['MCV1'])

        xval = np.arange(0, run_years, 1/12) + 1/24
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
            dose_dat[sidx, :] = np.array(data_brick[skey]['camp_cost'])
            mort_dat[sidx, :] = np.sum(age_dat[sidx, :, :]*mort_prob, axis=1)

        fidx = (pyr_mat[:, 0, 0] >= 0)

        pyr_mat_avg = np.mean(pyr_mat[fidx, :, :], axis=0)
        tpop_avg = np.sum(pyr_mat_avg, axis=1)
        tpop_xval = np.arange(len(tpop_avg))
        pops = np.interp(xyrs, tpop_xval, tpop_avg)
        mort_nrm = mort_dat[fidx, :]/pops*1e5
        dose_nrm = dose_dat[fidx, :]/pops*1e5

        # Figures
        fig01 = plt.figure(figsize=(8, 6))

        # Figures - Sims - Infections
        axs01 = fig01.add_subplot(1, 1, 1)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)

        axs01.set_ylabel('Monthly Mortality per-100k', fontsize=16)
        axs01.set_xlabel('MCV Coverage', fontsize=16)
        axs01.set_ylim(0, 4.0)
        axs01.set_xlim(0.1, 1.0)

        dval = np.round(np.mean(dose_nrm[fidx, -12:]))
        print('Avg annual SIA doses per-100k total pop: ', dval)

        yval = np.mean(mort_nrm[fidx, -12:], axis=1)/12.0
        print('Avg monthly mortality per-100k total pop: ', np.round(np.mean(yval), 3))

        print()

        axs01.plot(len(yval)*[mcv1_val], yval, '.', color='C0',
                   markeredgecolor=None, label=None)

        axs01.tick_params(axis='both', which='major', labelsize=14)

        plt.tight_layout()
        plt.savefig('fig_trends_{:s}01.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
