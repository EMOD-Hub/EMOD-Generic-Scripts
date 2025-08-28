# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))
from py_assets_common.emod_constants import NUM_SIMS, P_FILE, POP_PYR, EXP_V
from global_data import run_years, AGE_HIST_BINS

# *****************************************************************************

DIRNAMES = ['experiment_sweep_base',
            'experiment_sweep_base_MCV2',
            'experiment_sweep_base_SIAs',]

# *****************************************************************************


def make_fig():

    for dirname in DIRNAMES:

        # Sim outputs
        tpath = os.path.join('..', dirname)

        with open(os.path.join(tpath, 'data_brick.json')) as fid01:
            data_brick = json.load(fid01)

        with open(os.path.join(tpath, P_FILE)) as fid01:
            param_dict = json.load(fid01)

        nsims = int(param_dict[NUM_SIMS])

        inf_dat = np.zeros((nsims, 12*int(run_years)))
        age_dat = np.zeros((nsims, int(run_years), len(AGE_HIST_BINS)-1))
        pyr_mat = np.zeros((nsims, int(run_years)+1, 20))-1

        mcv1_vec = np.array(param_dict[EXP_V]['MCV1'])
        mcv1_lev = sorted(np.unique(mcv1_vec).tolist())

        num_charts = int(run_years//10)

        for skey in data_brick:
            if (not skey.isdigit()):
                continue

            sidx = int(skey)
            inf_dat[sidx, :] = np.array(data_brick[skey]['timeseries'])
            age_dat[sidx, :, :] = np.array(data_brick[skey]['age_data'])
            pyr_mat[sidx, :, :] = np.array(data_brick[skey][POP_PYR])

        fidx = (pyr_mat[:, 0, 0] >= 0)

        inf_yrs = np.zeros((nsims, int(run_years)))
        for k1 in range(int(run_years)):
            inf_yrs[:, k1] = np.sum(inf_dat[:, (k1*12):((k1+1)*12)], axis=1)

        bin_wid = np.diff(AGE_HIST_BINS)
        bin_loc = AGE_HIST_BINS[:-1] + bin_wid/2
        age_frac = age_dat/(inf_yrs[:, :, np.newaxis]+np.finfo(float).eps)
        age_frac_mean = np.mean(age_frac[fidx, :, :], axis=0)
        age_frac_std = np.std(age_frac[fidx, :, :], axis=0)

        # Figure
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
            axs01.set_xlim(min(AGE_HIST_BINS), max(AGE_HIST_BINS))
            axs01.text(8.3, 0.55, 'Year {:d}'.format(10*k1), fontsize=18)

        plt.tight_layout()
        plt.savefig('fig_agehist_{:s}_01.png'.format(dirname))
        plt.close()

        # Figure
        #fig01 = plt.figure(figsize=(8, 6))

        # x_vec = np.zeros(len(mcv1_lev))
        # y_vec = np.zeros(len(mcv1_lev))
        # age_dat_sum = np.sum(age_dat[:, -10:, :], axis=1)
        # inf_dat_sum = np.sum(inf_yrs[:, -10:], axis=1)

        # for k1 in range(len(mcv1_lev)):
            # gidx = fidx & (mcv1_vec==mcv1_lev[k1])
            # x_vec[k1] = mcv1_lev[k1]
            # y_vec[k1] = np.mean(np.sum(age_dat_sum[gidx, 6:9], axis=1)/inf_dat_sum[gidx])

        # axs01 = fig01.add_subplot(1, 1, 1)
        # plt.sca(axs01)

        # axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        # axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        # axs01.set_axisbelow(True)

        # axs01.plot(x_vec, y_vec)

        # axs01.set_ylabel('Burden Fraction: 6mo - 9mo', fontsize=16)
        # axs01.set_xlabel('MCV1 Coverage', fontsize=16)
        # axs01.set_xlim(0.20, 1.00)

        # plt.tight_layout()
        # plt.savefig('fig_ageburd_{:s}_01.png'.format(dirname))
        # plt.close()


    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
