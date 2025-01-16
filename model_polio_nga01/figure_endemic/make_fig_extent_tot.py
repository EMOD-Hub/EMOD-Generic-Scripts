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

from py_assets_common.emod_local_proc import shape_patch, shape_line

from global_data import base_year, init_ob_thresh

# *****************************************************************************

DIRNAMES = ['experiment_cVDPV2_NGA_100km_baseline_v01',
            'experiment_cVDPV2_NGA_100km_baseline_v02']

# *****************************************************************************


def make_fig():

    tpath = os.path.join('..', 'Assets', 'data','shapes_NGA00_COUNTRY.json')
    with open(tpath) as fid01:
        nga_shp00 = json.load(fid01)

    tpath = os.path.join('..', 'Assets', 'data','shapes_NGA02_LGA.json')
    with open(tpath) as fid01:
        nga_shp02 = json.load(fid01)

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
        run_years = int(param_dict[EXP_C]['run_years'])

        inf_data = np.zeros((n_sims, len(n_dict), t_vec.shape[0]))
        for sim_idx_str in data_brick:
            sim_idx = int(sim_idx_str)
            infmat = np.array(data_brick[sim_idx_str]['infmat'])
            inf_data[sim_idx, :, :] = infmat

        totinf = np.sum(inf_data, axis=1)
        cuminf = np.cumsum(totinf, axis=1)
        gidx = (cuminf[:, -1] >= init_ob_thresh)
        lgamat = (inf_data[gidx,:,:]>0)
        totlga = np.sum(lgamat, axis=1)

        #gidx = (cuminf[:, -1] > 1e6) & (cuminf[:, -1] < 2.5e6)
        #gidx = (totinf[:, -1] > 0) & (totlga[:, -240] > 15) & (totlga[:, -145] > 10)
        #print(np.argwhere(gidx))

        # Figure setup
        ax_pat = [run_years*[0], run_years*[0], run_years*[0],
                  list(range(1, run_years+1))]
        (fig01, axlist) = plt.subplot_mosaic(ax_pat, figsize=(2.5+1.5*run_years, 8))
        axs01 = axlist[0]

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)
        axs01.tick_params(axis='x', which='major', labelsize=18)
        axs01.tick_params(axis='y', which='major', labelsize=14)

        dvals = [0]+MO_DAYS*int(run_years)
        ticloc01 = np.cumsum(dvals) + t_vec[0]
        axs01.set_xticks(ticks=ticloc01, minor=True)

        ticloc02 = np.arange(0, int(run_years)+1) + t_vec[0]
        axs01.set_xticks(ticks=ticloc02)

        obp_lab = 'Fraction: {:5.3f}'.format(np.sum(gidx)/n_sims)
        axs01.text(0.05, 0.9, obp_lab, fontsize=14, transform = axs01.transAxes)

        yval2 = totlga
        yval1 = np.mean(yval2, axis=0)
        for k3 in range(yval2.shape[0]):
            #axs01.plot(t_vec, yval2[k3, :], '.', c='C0', alpha=0.1)
            axs01.plot(t_vec, yval2[k3, :], c='C0')
        #axs01.plot(t_vec, yval1, c='k', lw=3)

        axs01.set_ylabel('LGAs with Transmission', fontsize=14)
        axs01.set_xlim(t_vec[0], t_vec[-1])
        axs01.set_ylim(0, 250)

        nga0_prt = nga_shp00['AFRO:NIGERIA']['parts']
        nga0_pts = nga_shp00['AFRO:NIGERIA']['points']
        for k1 in range(run_years):
            axs01 = axlist[k1+1]
            axs01.axis('off')
            axs01.set_aspect('equal')
            shape_patch(axs01, nga0_pts, nga0_prt, clr=3*[0.9])
            yidx = (t_vec>=ticloc02[k1]) & (t_vec<ticloc02[k1+1])
            yrdat = np.max(lgamat[:, :, yidx], axis=2)
            yrdat = np.mean(yrdat, axis=0)
            for lga_name in n_dict:
                k2 = n_dict[lga_name]
                if(yrdat[k2] > 0):
                    nga2_prt = nga_shp02[lga_name]['parts']
                    nga2_pts = nga_shp02[lga_name]['points']
                    shape_patch(axs01, nga2_pts, nga2_prt,
                                clr=[1.0, 1.0-yrdat[k2], 1.0-yrdat[k2]])

        plt.tight_layout()
        plt.savefig('fig_extent_tot_{:s}_01.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
