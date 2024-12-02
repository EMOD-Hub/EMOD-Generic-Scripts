# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))
from py_assets_common.emod_constants import NUM_SIMS, P_FILE, D_FILE, \
                                            MO_DAYS, EXP_C
from py_assets_common.emod_local_proc import haversine_dist

from global_data import base_year
from refdat_location_admin02 import data_dict as ref_longlatref

# *****************************************************************************

DIRNAMES = ['experiment_cVDPV2_NGA_100km_baseline']

# *****************************************************************************


def make_fig():

    for dirname in DIRNAMES:

        # Sim outputs
        tpath = os.path.join('..', dirname)

        with open(os.path.join(tpath, D_FILE)) as fid01:
            data_brick = json.load(fid01)

        with open(os.path.join(tpath, P_FILE)) as fid01:
            param_dict = json.load(fid01)

        t_vec = np.array(data_brick.pop('t_vec'))/365 + base_year
        n_dict = data_brick.pop('node_names')
        n_sims = param_dict[NUM_SIMS]
        run_years = param_dict[EXP_C]['run_years']
        seed_inf_loc = param_dict[EXP_C]['seed_location']

        i_xy = ref_longlatref[seed_inf_loc]
        c_long = [ref_longlatref[nname][0] for nname in n_dict]
        c_lat = [ref_longlatref[nname][1] for nname in n_dict]
        dist_vec = haversine_dist(c_lat, c_long, i_xy[1], i_xy[0])

        inf_data = np.zeros((n_sims, len(n_dict), t_vec.shape[0]))
        for sim_idx_str in data_brick:
            sim_idx = int(sim_idx_str)
            infmat = np.array(data_brick[sim_idx_str]['infmat'])
            inf_data[sim_idx, :, :] = infmat

        totinf = np.sum(inf_data, axis=1)
        cuminf = np.cumsum(totinf, axis=1)
        dist_mat = dist_vec[np.newaxis, :, np.newaxis]
        lga_bool = (inf_data>0)*dist_mat
        max_dist = np.max(lga_bool, axis=1)

        gidx = (cuminf[:, -1] > 5000)
        gidx = (totinf[:, -1] > 0)

        # Figure setup
        fig01 = plt.figure(figsize=(8, 6))

        axs01 = fig01.add_subplot(1, 1, 1)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)
        axs01.tick_params(axis='x', which='major', labelsize=18)
        axs01.tick_params(axis='y', which='major', labelsize=14)

        dvals = [0]+MO_DAYS*int(run_years)
        ticloc = np.cumsum(dvals) + t_vec[0]
        axs01.set_xticks(ticks=ticloc, minor=True)

        ticloc = np.arange(0, int(run_years)+1) + t_vec[0]
        axs01.set_xticks(ticks=ticloc)

        obp_lab = 'Outbreak Probability: {:4.2f}'.format(np.sum(gidx)/n_sims)
        axs01.text(0.1, 0.9, obp_lab, fontsize=14, transform = axs01.transAxes)

        yval2 = max_dist[gidx, :]
        yval1 = np.mean(yval2, axis=0)
        for k3 in range(np.sum(gidx)):
            axs01.plot(t_vec, yval2[k3, :], '.', c='C0', alpha=0.1)
        axs01.plot(t_vec, yval1, c='k', lw=3)

        axs01.set_ylabel('Distance from Emergence (km)', fontsize=14)
        axs01.set_xlim(t_vec[0], t_vec[-1])
        axs01.set_ylim(0, 900)

        plt.tight_layout()
        plt.savefig('fig_extent_max_{:s}_01.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
