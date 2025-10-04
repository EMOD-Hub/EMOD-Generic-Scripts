# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))
from py_assets_common.emod_constants import EXP_V, CBR_VEC, \
                                            NUM_SIMS, P_FILE, POP_PYR, \
                                            POP_AGE_DAYS, R0_VEC, R0_TIME, \
                                            BASE_YEAR, D_FILE
from global_data import run_years, start_year, t_step_days, inf_prd_mean

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
        tref = (start_year-BASE_YEAR)*365
        ri_vec = np.array(param_dict[EXP_V]['RI_rate'])

        ri_lev = sorted(list(set(ri_vec.tolist())), reverse=True)
        pyr_mat = np.zeros((nsims, int(run_years)+1, 20))-1
        inf_mat = np.zeros((nsims, int(run_years), 20))
        r0_mat = np.zeros((nsims, int(run_years*365/t_step_days)))
        birth_mat = np.zeros((nsims, int(run_years)))

        for sim_idx_str in data_brick:
            sim_idx = int(sim_idx_str)
            sim_data = data_brick[sim_idx_str]
            pyr_mat[sim_idx, :, :] = np.array(sim_data[POP_PYR])
            inf_mat[sim_idx, :, :] = np.array(sim_data['inf_data'])
            birth_mat[sim_idx, :] = np.array(sim_data[CBR_VEC])

            r0_idx = (np.array(sim_data[R0_TIME])-tref)/t_step_days
            r0_idx = r0_idx.astype(int)
            r0_mat[sim_idx, r0_idx] = np.array(sim_data[R0_VEC])

        # Index for simulations with output
        fidx = (pyr_mat[:, 0, 0] >= 0)

        # Average population
        pyr_mat_avg = np.mean(pyr_mat[fidx, :, :], axis=0)
        pop_tot = np.sum(pyr_mat_avg, axis=1)
        pop_tot = np.diff(pop_tot)/2.0 + pop_tot[:-1]

        # Mean age-of-infection calculations
        pop_age_vec = np.array(POP_AGE_DAYS[:-1]) + np.diff(POP_AGE_DAYS)
        np_eps = np.finfo(float).eps
        tot_inf_mat = np.sum(inf_mat, axis=2) + np_eps
        avg_age = np.sum(inf_mat*pop_age_vec, axis=2)/tot_inf_mat

        # Figures
        fig01 = plt.figure(figsize=(24, 6))
        end_year = start_year+run_years
        dtstep01 = 1.0
        dtstep02 = t_step_days/365.0
        XDAT = np.arange(start_year, end_year, dtstep01) + dtstep01/2.0
        XDAT2 = np.arange(start_year, end_year, dtstep02) + dtstep02/2.0

        # Figures - Sims - Infections
        axs01 = fig01.add_subplot(1, 3, 1)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)

        ylab_str = 'Annual Rubella Infections per 100k Population'
        axs01.set_ylabel(ylab_str, fontsize=16)

        axs01.set_xlim(2020, 2050)
        axs01.set_ylim(0, 6000)

        axs01.set_yticks(ticks=np.arange(0, 6001, 1000))
        axs01.set_yticklabels(['0', '1k', '2k',
                               '3k', '4k', '5k', '6k'], fontsize=16)
        axs01.tick_params(axis='x', labelsize=16)

        for ri_val in ri_lev:
            ri_perc = 100*ri_val
            lab_str = 'RI = {:3d}%'.format(int(ri_perc))
            c_str = 'C{:d}'.format(int(ri_perc/20))
            gidx = (ri_vec == ri_val) & fidx
            inf_mat_avg = np.mean(inf_mat[gidx, :, :], axis=0)
            ydat = np.sum(inf_mat_avg, axis=1)/pop_tot*1e5

            axs01.plot(XDAT, ydat, label=lab_str, c=c_str)

        axs01.legend(fontsize=14)

        # Figures - Sims - Average Age of Infection
        axs01 = fig01.add_subplot(1, 3, 2)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)

        axs01.set_ylabel('Mean Age at Infection (yrs)', fontsize=16)

        axs01.set_xlim(2020, 2050)
        axs01.set_ylim(0.0, 35.0)

        axs01.tick_params(axis='x', labelsize=16)
        axs01.tick_params(axis='y', labelsize=16)

        for ri_val in ri_lev:
            ri_perc = 100*ri_val
            lab_str = 'RI = {:3d}%'.format(int(ri_perc))
            c_str = 'C{:d}'.format(int(ri_perc/20))
            gidx = (ri_vec == ri_val) & fidx
            ydat = 0.0*XDAT
            for k1 in range(ydat.shape[0]):
                nzgidx = (avg_age[:, k1] > 0) & gidx
                ydat[k1] = np.mean(avg_age[nzgidx, k1], axis=0)/365.0

            axs01.plot(XDAT, ydat, label=lab_str, c=c_str)

        # Figures - Sims - Estimated R0
        axs01 = fig01.add_subplot(1, 3, 3)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)

        axs01.set_ylabel('Estimated R$_{0}$ Value', fontsize=16)

        axs01.set_xlim(2020, 2050)
        axs01.set_ylim(0.0, 10.0)

        axs01.tick_params(axis='x', labelsize=16)
        axs01.tick_params(axis='y', labelsize=16)

        for ri_val in ri_lev:
            ri_perc = 100*ri_val
            lab_str = 'RI = {:3d}%'.format(int(ri_perc))
            c_str = 'C{:d}'.format(int(ri_perc/20))
            gidx = (ri_vec == ri_val) & fidx
            ydat = 0.0*XDAT2
            for k1 in range(ydat.shape[0]):
                nzgidx = (r0_mat[:, k1] > 0) & gidx
                if (np.sum(nzgidx)):
                    ydat[k1] = np.mean(r0_mat[nzgidx, k1], axis=0)*inf_prd_mean

            axs01.plot(XDAT2, ydat, label=lab_str, c=c_str)

        plt.tight_layout()
        plt.savefig('fig_aai_{:s}01.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
