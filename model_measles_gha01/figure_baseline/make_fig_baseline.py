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
from py_assets_common.emod_constants import NUM_SIMS, P_FILE, D_FILE, \
                                            EXP_C

from global_data import base_year

# *****************************************************************************

DIRNAMES = ['experiment_gha_base01']

# *****************************************************************************


def make_fig():

    YMAX = 350
    CCUT = -50e3

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
            sim_obj = dbrick[sim_idx_str]

            inf_data[sim_idx, :] = np.array(sim_obj['timeseries'])
            scale_vec[sim_idx, 0] = sim_obj['rep_rate']
            cal_vec[sim_idx] = sim_obj['cal_val']

        gidx = (scale_vec[:,0] >= 0)
        gidx = gidx & (cal_vec > CCUT)

        # Figure
        fig01 = plt.figure(figsize=(24, 6))

        # Timeseries
        axs01 = fig01.add_subplot(1, 3, 1)
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

        xval = np.array(tvals)/365+base_year

        inf_data_sort = np.sort(inf_data[gidx, :]*scale_vec[gidx], axis=0)
        for patwid in [0.475, 0.375, 0.25]:
            xydat = np.zeros((2*inf_data_sort.shape[1], 2))
            xydat[:, 0] = np.hstack((xval, xval[::-1]))
            tidx = int((0.5-patwid)*inf_data_sort.shape[0])
            xydat[:,1] = np.hstack((inf_data_sort[tidx, :],
                                    inf_data_sort[-tidx, ::-1]))

            poly_shp = patch.Polygon(xydat, facecolor='C0',
                                     alpha=0.7-patwid, edgecolor=None)
            axs01.add_patch(poly_shp)

        yval = np.mean(inf_data[gidx]*scale_vec[gidx], axis=0)
        axs01.plot(xval, yval, color='C0', linewidth=2)

        axs01.set_xlim(2010, 2020)
        axs01.set_ylim(0, YMAX)

        tpath = os.path.join('..', 'Assets', 'data', 'GHA_MCV_SIA.json')
        with open(tpath) as fid01:
            ref_dat_sia = json.load(fid01)

        for sia_name in ref_dat_sia:
            xval = ref_dat_sia[sia_name]['date']
            sia01 = np.array([xval, xval])/365+base_year
            lstyle = ':'
            lwidth = 3
            llabel = 'SIA'
            if (ref_dat_sia[sia_name]['age_yr_max'] > 5.0):
                lstyle = '--'
                lwidth = 5
                llabel = 'RCV Catch-up'
            axs01.plot(sia01, [0, YMAX], color='c',
                       linewidth=lwidth, ls=lstyle)
            axs01.text(sia01[0]+0.2, 0.9*YMAX, llabel, fontsize=13)

        # Calibration score
        axs01 = fig01.add_subplot(1, 3, 2)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)

        axs01.hist(cal_vec, edgecolor='k', bins=np.arange(-8e3, -2e3, 200))

        # Reporting rate
        axs01 = fig01.add_subplot(1, 3, 3)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)

        axs01.hist(1/scale_vec, edgecolor='k', bins=100)

        plt.tight_layout()
        plt.savefig('fig_baseline_01.png')
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
