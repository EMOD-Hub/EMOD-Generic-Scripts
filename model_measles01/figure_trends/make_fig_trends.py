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

        nsims = int(param_dict[NUM_SIMS])

        inf_dat = np.zeros((nsims, 12*int(run_years)))
        age_dat = np.zeros((nsims, int(run_years), len(AGE_HIST_BINS)-1))
        mort_dat = np.zeros((nsims, int(run_years)))
        pyr_mat = np.zeros((nsims, int(run_years)+1, 20))-1

        mcv1_vec = np.array(param_dict[EXP_V]['MCV1'])
        mcv1_lev = sorted(np.unique(mcv1_vec).tolist())

        mcv1_age_vec = np.array(param_dict[EXP_V]['MCV1_age'])
        mcv1_age_lvl = np.unique(mcv1_age_vec).tolist()

        mcv2_lvl = param_dict[EXP_C]['MCV2']

        sia_year = param_dict[EXP_C]['sia_start_year']
        if ('sia_min_age_yr' in param_dict[EXP_V]):
            sia_min_age_yr_vec = np.array(param_dict[EXP_V]['sia_min_age_yr'])
        else:
            sia_min_age_yr_val = param_dict[EXP_C]['sia_min_age_yr']
            sia_min_age_yr_vec = np.array(nsims*[sia_min_age_yr_val])
        sia_min_age_yr_lvl = np.unique(sia_min_age_yr_vec).tolist()

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
            mort_dat[sidx, :] = np.sum(age_dat[sidx, :, :]*mort_prob, axis=1)

        fidx = (pyr_mat[:, 0, 0] >= 0)

        pyr_mat_avg = np.mean(pyr_mat[fidx, :, :], axis=0)
        tpop_avg = np.sum(pyr_mat_avg, axis=1)
        tpop_xval = np.arange(len(tpop_avg))
        pops = np.interp(xyrs, tpop_xval, tpop_avg)
        mort_nrm = mort_dat[fidx, :]/pops*1e5

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
        axs01.set_xlim(0.2, 1.0)

        for k1 in range(len(mcv1_age_lvl)):
            mcv1_age_val = mcv1_age_lvl[k1]
            idx02 = (mcv1_age_vec == mcv1_age_val)
            for k2 in range(len(sia_min_age_yr_lvl)):
                sia_min_age_yr_val = sia_min_age_yr_lvl[k2]
                idx03 = (sia_min_age_yr_vec == sia_min_age_yr_val)

                lstyle = '-'
                if (sia_min_age_yr_val < 0.7):
                    axs01.set_ylim(0, 1.8)
                    lstyle = '--'

                tidx = (fidx & idx02 & idx03)
                xval = mcv1_vec[tidx]
                xval2 = np.arange(0.2, 1.001, 0.01)
                yval = np.mean(mort_nrm[tidx, -10:], axis=1)/12.0
                pcoef = np.polyfit(xval, yval, 5)
                yval2 = np.polyval(pcoef, xval2)

                cval = 'C{:d}'.format(k1)

                mcv1_mo = int(np.round(mcv1_age_val/365*12))
                lstr = 'MCV1 {:>2d}mo'.format(mcv1_mo)
                if (mcv2_lvl > 0):
                    lstr = lstr + '+ MCV2 15mo'
                if (sia_year < run_years):
                    sia_mo = int(sia_min_age_yr_val*12)
                    lstr = lstr + '+ SIAs >{:d}mo'.format(sia_mo)

                axs01.plot(xval, yval, '.', alpha=0.1, color=cval,
                           markeredgecolor=None, label=None)
                axs01.plot(xval2, yval2, lstyle, alpha=1.0, color=cval,
                           label=lstr)

        axs01.tick_params(axis='both', which='major', labelsize=14)
        axs01.legend(fontsize=16)

        plt.tight_layout()
        plt.savefig('fig_trends_{:s}01.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
