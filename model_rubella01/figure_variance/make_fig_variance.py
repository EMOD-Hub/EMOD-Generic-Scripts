# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))
from py_assets_common.emod_local_proc import crs_proc
from py_assets_common.emod_constants import EXP_C, EXP_V, CBR_VEC, \
                                            NUM_SIMS, P_FILE, D_FILE, POP_PYR
from global_data import run_years, start_year

# *****************************************************************************

DIRNAMES = ['experiment_SSA_sweepRI_popEQL_noSIAs',
            'experiment_SSA_sweepRI_popMED_noSIAs']

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
        ss_demog = param_dict[EXP_C]['steady_state_demog']
        demog_set = param_dict[EXP_C]['demog_set']
        ri_vec = np.array(param_dict[EXP_V]['RI_rate'])

        ri_lev = sorted(list(set(ri_vec.tolist())))
        pyr_mat = np.zeros((nsims, int(run_years)+1, 20))-1
        inf_mat = np.zeros((nsims, int(run_years), 20))
        birth_mat = np.zeros((nsims, int(run_years)))

        for sim_idx_str in data_brick:
            sim_idx = int(sim_idx_str)
            sim_data = data_brick[sim_idx_str]
            pyr_mat[sim_idx, :, :] = np.array(sim_data[POP_PYR])
            inf_mat[sim_idx, :, :] = np.array(sim_data['inf_data'])
            birth_mat[sim_idx, :] = np.array(sim_data[CBR_VEC])

        # Index for simulations with output
        fidx = (pyr_mat[:, 0, 0] >= 0)

        # Average population
        pyr_mat_avg = np.mean(pyr_mat[fidx, :, :], axis=0)
        pop_tot = np.sum(pyr_mat_avg, axis=1)
        pop_tot = np.diff(pop_tot)/2.0 + pop_tot[:-1]

        # CRS calculations
        XDAT = np.arange(start_year, start_year+run_years) + 0.5
        fname = 'fert_dat_{:s}.csv'.format(demog_set)
        fnabs = os.path.abspath(os.path.join('..', 'Assets', 'data', fname))
        (frt_brth, crs_prob_vec) = crs_proc(fnabs, XDAT, pyr_mat_avg, ss_demog)

        # Normalize timeseries required for CRS calculation
        brth_vec = np.mean(birth_mat[fidx, :], axis=0)
        norm_crs_timevec = brth_vec/frt_brth

        # Figures
        fig01 = plt.figure(figsize=(16, 6))

        # Figures - Sims - CRS
        axs01 = fig01.add_subplot(1, 2, 1)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)

        axs01.set_ylabel('Annual CRS Burden per 1k Births', fontsize=16)

        axs01.set_xlim(2020, 2060)
        axs01.set_ylim(0.0, 5.0)

        axs01.set_yticks(ticks=np.arange(0, 5.1, 0.5))
        axs01.set_yticklabels(['0', '', '1', '', '2', '',
                               '3', '', '4', '', '5'], fontsize=16)
        axs01.tick_params(axis='x', labelsize=16)

        for k1 in range(len(ri_lev)):
            ri_val = ri_lev[k1]
            if (ri_val not in [0.0, 0.6]):
                continue
            gidx = (ri_vec == ri_val) & fidx
            inf_mat_avg = np.mean(inf_mat[gidx, :, :], axis=0)
            crs_mat = inf_mat_avg*np.transpose(crs_prob_vec)
            ydat = np.sum(crs_mat, axis=1)/brth_vec*norm_crs_timevec*1e3

            label_str = 'RI = {:3d}%'.format(int(100*ri_val))
            color_str = 'C{:d}'.format(k1)
            axs01.plot(XDAT, ydat, c=color_str, label=label_str)

        # Figures - Sims - CRS Histogram
        axs01 = fig01.add_subplot(1, 2, 2)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)

        xlab_str = 'Annual Rubella Burden per 1k Births: 2050 to 2060'
        axs01.set_xlabel(xlab_str, fontsize=16)
        axs01.set_ylabel('Probability', fontsize=16)

        xmaxv = 5
        xbinv = 0.1
        axs01.set_xlim(0.0, xmaxv)
        axs01.set_ylim(0.0, 4.0)

        axs01.tick_params(axis='x', labelsize=16)
        axs01.tick_params(axis='y', labelsize=16)

        for k1 in range(len(ri_lev)):
            ri_val = ri_lev[k1]
            if (ri_val not in [0.0, 0.6]):
                continue
            gidx = (ri_vec == ri_val) & fidx
            inf_mat_sub = inf_mat[gidx, :, :]
            crs_mat = inf_mat_sub*np.transpose(crs_prob_vec)
            ydat = np.sum(crs_mat, axis=2)/brth_vec*norm_crs_timevec*1e3
            ydat = ydat[:, -10:].flatten()

            label_str = 'RI = {:3d}%'.format(int(100*ri_val))
            color_str = 'C{:d}'.format(k1)
            axs01.hist(ydat, bins=np.arange(0, xmaxv+xbinv, xbinv),
                       density=True, alpha=0.7, color=color_str,
                       label=label_str)
            print(ri_val, np.mean(ydat), np.percentile(ydat, [5, 50, 95]))

        axs01.legend(fontsize=14)

        plt.tight_layout()
        plt.savefig('fig_var_{:s}01.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
