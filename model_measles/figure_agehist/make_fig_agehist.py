# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))

from py_assets_common.emod_constants import NUM_SIMS, P_FILE, POP_PYR, EXP_V, \
                                            D_FILE
from global_data import run_years, AGE_HIST_BINS, IHME_MORT_X, IHME_MORT_Y

# *****************************************************************************

DIRNAMES = ['experiment_sweep_base']

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

        age_dat = np.zeros((NSIMS, int(run_years), len(AGE_HIST_BINS)-1))
        mort_dat = np.zeros((NSIMS, int(run_years), len(AGE_HIST_BINS)-1))
        pyr_mat = np.zeros((NSIMS, int(run_years)+1, 20))-1

        mcv1_vec = np.array(param_dict[EXP_V]['MCV1'])
        mcv1_lev = sorted(np.unique(mcv1_vec).tolist())

        bin_wid = np.diff(AGE_HIST_BINS)
        bin_loc = AGE_HIST_BINS[:-1] + bin_wid/2
        mort_prob = np.interp(bin_loc, IHME_MORT_X, IHME_MORT_Y)

        for skey in data_brick:
            if (not skey.isdigit()):
                continue

            sidx = int(skey)
            age_dat[sidx, :, :] = np.array(data_brick[skey]['age_data'])
            mort_dat[sidx, :, :] = age_dat[sidx, :, :]*mort_prob
            pyr_mat[sidx, :, :] = np.array(data_brick[skey][POP_PYR])

        fidx = (pyr_mat[:, 0, 0] >= 0)

        inf_yrs = np.sum(age_dat, axis=2)
        age_frac = age_dat/(inf_yrs[:, :, np.newaxis]+np.finfo(float).eps)
        age_frac_mean = np.mean(age_frac[fidx, :, :], axis=0)
        age_frac_std = np.std(age_frac[fidx, :, :], axis=0)

        mort_yrs = np.sum(mort_dat, axis=2)
        mort_frac = mort_dat/(mort_yrs[:, :, np.newaxis]+np.finfo(float).eps)
        mort_frac_mean = np.mean(mort_frac[fidx, :, :], axis=0)
        mort_frac_std = np.std(mort_frac[fidx, :, :], axis=0)

        # Figure
        num_charts = int(run_years//10)
        fig01 = plt.figure(figsize=(8*num_charts, 6))

        for k1 in range(1, num_charts+1):

            axs01 = fig01.add_subplot(1, num_charts, k1)
            plt.sca(axs01)

            axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
            axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
            axs01.set_axisbelow(True)

            bin_hgt = age_frac_mean[10*k1-1, :]/bin_wid
            bin_std = age_frac_std[10*k1-1, :]/bin_wid
            axs01.bar(bin_loc, bin_hgt, width=bin_wid, edgecolor='k',
                      yerr=bin_std)

            axs01.set_ylabel('Fraction', fontsize=16)
            axs01.set_xlabel('Age at Infection (yrs)', fontsize=16)
            axs01.set_ylim(0.0, 0.65)
            axs01.set_xlim(min(AGE_HIST_BINS), 10.0)
            axs01.text(8.3, 0.55, 'Year {:d}'.format(10*k1), fontsize=18)

        plt.tight_layout()
        plt.savefig('fig_agehist_inf_{:s}01.png'.format(dirname))
        plt.close()

        # Figure
        num_charts = int(run_years//10)
        fig01 = plt.figure(figsize=(8*num_charts, 6))

        for k1 in range(1, num_charts+1):

            axs01 = fig01.add_subplot(1, num_charts, k1)
            plt.sca(axs01)

            axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
            axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
            axs01.set_axisbelow(True)

            bin_hgt = mort_frac_mean[10*k1-1, :]/bin_wid
            bin_std = mort_frac_std[10*k1-1, :]/bin_wid
            axs01.bar(bin_loc, bin_hgt, width=bin_wid, edgecolor='k',
                      yerr=bin_std)

            axs01.set_ylabel('Fraction', fontsize=16)
            axs01.set_xlabel('Age at Infection (yrs)', fontsize=16)
            axs01.set_ylim(0.0, 0.65)
            axs01.set_xlim(min(AGE_HIST_BINS), 10.0)
            axs01.text(8.3, 0.55, 'Year {:d}'.format(10*k1), fontsize=18)

        plt.tight_layout()
        plt.savefig('fig_agehist_mort_{:s}01.png'.format(dirname))
        plt.close()


    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
