# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))
from py_assets_common.emod_constants import NUM_SIMS, P_FILE, D_FILE, EXP_C

from py_assets_common.emod_local_proc import shape_patch, shape_line

from global_data import base_year, init_ob_thresh

# *****************************************************************************

DIRNAMES = ['experiment_cVDPV2_NGA_100km_baseline']

# *****************************************************************************


def make_fig():

    dy_init = 1
    dy_end = 0

    tpath = os.path.join('..', 'Assets', 'data','shapes_NGA00_COUNTRY.json')
    with open(tpath) as fid01:
        nga_shp00 = json.load(fid01)

    tpath = os.path.join('..', 'Assets', 'data','shapes_NGA02_LGA.json')
    with open(tpath) as fid01:
        nga_shp02 = json.load(fid01)

    tpath = os.path.join('..', 'Assets', 'data','NGA_epi.csv')
    ref_dat_mo = np.loadtxt(tpath, dtype=int, delimiter=',')
    tvec_ref = np.arange(2016, 2026, 1/12) + 1/24

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
        year_init = int(param_dict[EXP_C]['start_year']) + dy_init
        run_years = int(param_dict[EXP_C]['run_years']) - dy_init - dy_end
        tbool = (t_vec >= year_init) & (t_vec < (year_init+run_years))
        tbool_ref = (tvec_ref >= year_init) & (tvec_ref < (year_init+run_years))

        inf_data = np.zeros((n_sims, len(n_dict), t_vec.shape[0]))
        for sim_idx_str in data_brick:
            sim_idx = int(sim_idx_str)
            infmat = np.array(data_brick[sim_idx_str]['infmat'])
            inf_data[sim_idx, :, :] = infmat

        totinf = np.sum(inf_data, axis=1)
        cuminf = np.cumsum(totinf, axis=1)
        gidx = (cuminf[:, -1] >= init_ob_thresh)

        #cumlga = np.cumsum(inf_data, axis=2)[:,:,-1]
        #gidx = gidx & (np.sum(cumlga[:, 141:167], axis=1)>0)
        #gidx = gidx & (np.max(totinf[:,:150], axis=1) < 50000)
        #gidx = gidx & (np.max(totinf[:,:150], axis=1) > 3000)
        #gidx = gidx & (totinf[:, -1] > 0) #& (cuminf[:, -1] > 900e3) & (cuminf[:, -1] < 1000e3)
        #gidx = gidx & (cuminf[:, -75] < 100e3) #& (cuminf[:, -74] > 40e3) & (cuminf[:, -1] < 60e3)
        #gidx = gidx & (cuminf[:, -50] < 120e3) #& (cuminf[:, -74] > 40e3) & (cuminf[:, -1] < 60e3)
        #gidx = gidx & (cuminf[:, -1] > 400e3)
        #gidx = gidx & ((cuminf[:, -1] - cuminf[:, -73]) > 45e3)
        #gidx = gidx & ((cuminf[:, -50] - cuminf[:, -73]) > 17e3)
        #print(np.argwhere(gidx))
        #gidx = gidx & (np.array(list(range(n_sims))) >= 200) & (np.array(list(range(n_sims))) < 500)
        gidx = (np.array(list(range(n_sims))) == 39)
        #print(cuminf[gidx, -1])

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

        ticloc01 = np.arange(0, int(run_years)+0.001) + t_vec[tbool][0]
        axs01.set_xticks(ticks=ticloc01)

        #obp_lab = 'Fraction: {:5.3f}'.format(np.sum(gidx)/n_sims)
        #axs01.text(0.05, 0.9, obp_lab, fontsize=14, transform = axs01.transAxes)

        yval1 = totinf[gidx]/1000
        yval2 = np.mean(yval1, axis=0)
        for k3 in range(yval1.shape[0]):
            axs01.plot(t_vec[tbool], yval1[k3, tbool], '.', c='C0', alpha=0.1)
            #axs01.plot(t_vec[tbool], yval1[k3, tbool])#, c='C0')
        axs01.plot(t_vec[tbool], yval2[tbool], c='k', lw=3)

        axs01.set_ylabel('Simulated Incidence (thousands)', fontsize=18)
        axs01.set_xlim(t_vec[tbool][0], t_vec[tbool][-1]+0.02)
        axs01.set_yscale('symlog', linthresh=25)
        axs01.set_ylim(0, 100)
        ticloc02 = list(range(0,31,5))+list(range(40,101,10))
        ticlab02 = [str(val) for val in ticloc02]
        ticlab02 = ticlab02[:9]+4*['']+[ticlab02[-1]]
        axs01.set_yticks(ticks=ticloc02)
        axs01.set_yticklabels(ticlab02)

        axs02 = axs01.twinx()
        axs02.bar(tvec_ref[tbool_ref], ref_dat_mo[tbool_ref], width=1/12,
                  alpha=0.2, facecolor='C0', edgecolor=None)
        axs02.set_xlim(t_vec[tbool][0], t_vec[tbool][-1]+0.02)
        axs02.set_yscale('symlog', linthresh=100)
        axs02.set_ylim(0, 400)
        ticloc02 = list(range(0,121,20))+list(range(160,401,40))
        ticlab02 = [str(val) for val in ticloc02]
        ticlab02 = ticlab02[:9]+4*['']+[ticlab02[-1]]
        axs02.set_yticks(ticks=ticloc02)
        axs02.set_yticklabels(ticlab02)
        axs02.tick_params(axis='y', which='major', labelsize=14)
        axs02.set_ylabel('Monthly AFP Cases', fontsize=18)

        nga0_prt = nga_shp00['AFRO:NIGERIA']['parts']
        nga0_pts = nga_shp00['AFRO:NIGERIA']['points']
        lgamat = (inf_data[gidx,:,:]>0)
        for k1 in range(run_years):
            axs01 = axlist[k1+1]
            axs01.axis('off')
            axs01.set_aspect('equal')
            shape_patch(axs01, nga0_pts, nga0_prt, clr=3*[0.9])
            yidx = (t_vec>=ticloc01[k1]) & (t_vec<ticloc01[k1+1])
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
        plt.savefig('fig_extent_working_{:s}_01.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
