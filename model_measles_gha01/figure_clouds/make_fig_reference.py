# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import matplotlib as mpl

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))
from py_assets_common.emod_analysis import norpois_opt
from py_assets_common.emod_constants import NUM_SIMS, P_FILE, D_FILE, \
                                            EXP_C

# *****************************************************************************

DIRNAMES = ['experiment_gha_base01']

# *****************************************************************************


def make_fig():

    YMAX = 350

    for dirname in DIRNAMES:

        # Sim outputs
        tpath = os.path.join('..', dirname)

        with open(os.path.join(tpath, P_FILE)) as fid01:
            param_dict = json.load(fid01)

        with open(os.path.join(tpath, D_FILE)) as fid01:
            dbrick = json.load(fid01)

        nsims = int(param_dict[NUM_SIMS])
        tvals = dbrick.pop('tstamps')

        inf_data = np.zeros((nsims, len(tvals)))
        scale_vec = np.zeros((nsims, 1))-1
        cal_vec = np.zeros(nsims)

        for sim_idx_str in dbrick:
            if (not sim_idx_str.isdigit()):
                continue
            sim_idx = int(sim_idx_str)

            inf_data[sim_idx, :] = np.array(dbrick[sim_idx_str]['timeseries'])
            #scale_vec[sim_idx, 0] = dbrick[sim_idx_str]['rep_rate']
            #cal_vec[sim_idx] = dbrick[sim_idx_str]['cal_val']

        fidx = (scale_vec[:,0] >= 0)

        # Figure
        fig01 = plt.figure(figsize=(8, 6))

        axs01 = fig01.add_subplot(1, 1, 1)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)

        tpath = os.path.join('..', 'Assets', 'data', 'GHA_epi.json')
        with open(tpath) as fid01:
            ref_dat = json.load(fid01)

        ref_cases = np.array(ref_dat['cases_mo'])
        ref_start = ref_dat['start_year']
        ref_years = int(np.ceil(ref_cases.shape[0]/12))
        ref_dates = np.arange(ref_start, ref_start+ref_years, 1/12) + 1/24

        axs01.bar(ref_dates, ref_cases, color='r', linewidth=0.5,
                  edgecolor='k', width=0.9/12)
        axs01.set_ylabel('Monthly Cases', fontsize=16)

        yval = np.mean(inf_data, axis=0)
        xval = np.array(tvals)/365+1900
        axs01.plot(xval, yval, color='C0', linewidth=2)

        axs01.set_xlim(2010, 2020)
        axs01.set_ylim(0, YMAX)

        ticloc = [2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020]
        ticlab = ['', '', '', '', '', '', '', '', '', '', '']
        axs01.set_xticks(ticks=ticloc)
        axs01.set_xticklabels(ticlab)
        for k1 in ticloc[:-1]:
            axs01.text(k1+0.5, -0.04*YMAX, str(k1), fontsize=11, ha='center')

        sia01 = np.array([40460, 40460])/365+1900
        axs01.plot(sia01, [0, YMAX], color='c', linewidth=3, ls=':')
        axs01.text(sia01[0]+0.2, 0.9*YMAX, 'SIA', fontsize=13)

        sia02 = np.array([41508, 41508])/365+1900
        axs01.plot(sia02, [0, YMAX], color='c', linewidth=5, ls='--')
        axs01.text(sia02[0]+0.2, 0.9*YMAX, 'RCV Catch-up', fontsize=13)

        sia03 = np.array([43365, 43365])/365+1900
        axs01.plot(sia03, [0, YMAX], color='c', linewidth=3, ls=':')
        axs01.text(sia03[0]+0.2, 0.9*YMAX, 'SIA', fontsize=13)

        plt.tight_layout()
        plt.savefig('fig_reference_01.png')
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
