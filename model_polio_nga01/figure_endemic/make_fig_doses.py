# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))
from py_assets_common.emod_constants import NUM_SIMS, P_FILE, D_FILE, EXP_C, \
                                            CAMP_COST

from py_assets_common.emod_local_proc import shape_patch, shape_line

from global_data import base_year, init_ob_thresh, targ_adm00

# *****************************************************************************

DIRNAMES = [
            ('experiment_cVDPV2_NGA_100km_baseline', 0),
            #('experiment_cVDPV2_NGA_100km_baseline_ob01', 0),
            #('experiment_cVDPV2_NGA_100km_baseline_ob02', 0),
            ]

# *****************************************************************************


def make_fig():

    dy_init = 0
    dy_end = 0

    tpath = os.path.join('..', 'Assets', 'data','shapes_adm00.json')
    with open(tpath) as fid01:
        adm00_shp = json.load(fid01)

    tpath = os.path.join('..', 'Assets', 'data','shapes_adm02.json')
    with open(tpath) as fid01:
        adm02_shp = json.load(fid01)

    for dir_tup in DIRNAMES:
 
        dirname = dir_tup[0]
        fig_clr = 'C{:d}'.format(dir_tup[1])

        # Sim outputs
        tpath = os.path.join('..', dirname)

        with open(os.path.join(tpath, 'data_brick.json')) as fid01:
            data_brick = json.load(fid01)

        with open(os.path.join(tpath, P_FILE)) as fid01:
            param_dict = json.load(fid01)

        t_vec = np.array(data_brick['t_vec'])/365 + base_year
        n_dict = data_brick['node_names']
        n_sims = param_dict[NUM_SIMS]
        year_init = int(param_dict[EXP_C]['start_year']) + dy_init
        run_years = int(param_dict[EXP_C]['run_years']) - dy_init - dy_end
        tbool = (t_vec >= year_init) & (t_vec < (year_init+run_years))

        inf_data = np.zeros((n_sims, len(n_dict), t_vec.shape[0]))
        sia_data = np.zeros((n_sims, t_vec.shape[0]))
        for sim_idx_str in data_brick:
            if (not sim_idx_str.isdecimal()):
                continue
            sim_idx = int(sim_idx_str)
            infmat = np.array(data_brick[sim_idx_str]['infmat'])
            inf_data[sim_idx, :, :] = infmat
            siavec = np.array(data_brick[sim_idx_str][CAMP_COST])
            sia_data[sim_idx, :] = siavec

        totinf = np.sum(inf_data, axis=1)
        cuminf = np.cumsum(totinf, axis=1)
        gidx = (cuminf[:, -1] >= init_ob_thresh)
        #gidx = gidx & (cuminf[:, -1] > 900e3) #& (cuminf[:, -1] < 180e3)
        #gidx = gidx & (cuminf[:, -1] > 150e3)
        #gidx = gidx & (np.max(totinf, axis=1) < 6e3)
        #gidx = gidx & (totinf[:, -1] > 0) #& (cuminf[:, -1] < 180e3) #& (cuminf[:, -1] > 120e3)
        #gidx = gidx & (np.array(list(range(n_sims))) == 55) #& (np.array(list(range(n_sims))) < 900) #104
        #gidx = gidx & (np.array(list(range(n_sims))) > 225) #& (np.array(list(range(n_sims))) <= 225)

        print(np.sum(gidx))
        #print(np.argwhere(gidx))
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

        yval1 = sia_data[gidx]/1e6
        yval2 = np.mean(yval1, axis=0)
        for k3 in range(yval1.shape[0]):
            axs01.plot(t_vec[tbool], yval1[k3, tbool], '.', c=fig_clr, alpha=0.1)
        axs01.plot(t_vec[tbool], yval2[tbool], c='k', lw=3)

        axs01.set_ylabel('Cumulative SIA Vaccine Doses (millions)', fontsize=16)
        axs01.set_xlim(t_vec[tbool][0], t_vec[tbool][-1]+0.02)
        axs01.set_ylim(0, 120)
        #ticloc02 = list(range(0,31,5))+list(range(40,101,10))
        #ticlab02 = [str(val) for val in ticloc02]
        #ticlab02 = ticlab02[:9]+4*['']+[ticlab02[-1]]
        #axs01.set_yticks(ticks=ticloc02)
        #axs01.set_yticklabels(ticlab02)

        (xmin, xmax, ymin, ymax) = (180, -180, 90, -90)
        for cname in targ_adm00:
            pts = np.array(adm00_shp[cname]['points'])
            xmin = np.floor(min(xmin, np.min(pts[:,0])))
            xmax = np.ceil(max(xmax, np.max(pts[:,0])))
            ymin = np.floor(min(ymin, np.min(pts[:,1])))
            ymax = np.ceil(max(ymax, np.max(pts[:,1])))
        xydelt = max(xmax-xmin, ymax-ymin)

        #adm02_mat02 = (inf_data[gidx,:,:]>1)
        for k1 in range(run_years):
            axs01 = axlist[k1+1]
            axs01.axis('off')
            axs01.set_xlim(xmin, xmin+xydelt)
            axs01.set_ylim(ymin, xmin+xydelt)
        #    for cname in targ_adm00:
        #        adm00_prt = adm00_shp[cname]['parts']
        #        adm00_pts = adm00_shp[cname]['points']
        #        shape_patch(axs01, adm00_pts, adm00_prt, clr=3*[0.9])
        #        shape_line(axs01, adm00_pts, adm00_prt, wid=0.2)
        #    yidx = (t_vec>=ticloc01[k1]) & (t_vec<ticloc01[k1+1])
        #    yrdat02 = np.max(adm02_mat02[:, :, yidx], axis=2)
        #    yrdat02 = np.mean(yrdat02, axis=0)
        #    for adm02_name in n_dict:
        #        k2 = n_dict[adm02_name]
        #        if (yrdat02[k2] > 0):
        #            adm02_prt = adm02_shp[adm02_name]['parts']
        #            adm02_pts = adm02_shp[adm02_name]['points']
        #            clr_val = [1.0, 1.0, 1.0]
        #            shape_patch(axs01, adm02_pts, adm02_prt, clr=clr_val)

        plt.tight_layout()
        plt.savefig('fig_doses_{:s}_02.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
