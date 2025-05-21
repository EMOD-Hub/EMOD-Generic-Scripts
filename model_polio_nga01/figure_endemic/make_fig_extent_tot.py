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
                                            MO_DAYS

from py_assets_common.emod_analysis import norpois_opt, norpois_vec
from py_assets_common.emod_local_proc import shape_patch, shape_line

from global_data import base_year, init_ob_thresh, targ_adm00, t_step_days

# *****************************************************************************

DIRNAMES = [
            ('experiment_cVDPV2_NGA_100km_baseline', 0),
            #('experiment_cVDPV2_NGA_100km_baseline__02', 0),
            #('experiment_cVDPV2_NGA_100km_baseline_ob02', 0),
            #('experiment_cVDPV2_NGA_100km_baseline_RI', 4),
            #('experiment_cVDPV2_NGA_100km_baseline_SIA01', 1),
            #('experiment_cVDPV2_NGA_100km_baseline_SIA01N', 7),
            #('experiment_cVDPV2_NGA_100km_baseline_SIA02', 7),
            #('experiment_cVDPV2_NGA_100km_baseline_SIA02N', 7),
            #('experiment_cVDPV2_NGA_100km_baseline_RI_SIA01', 2),
            #('experiment_cVDPV2_NGA_100km_baseline_RI_SIA01N', 5),
            #('experiment_cVDPV2_NGA_100km_baseline_RI_SIA02', 8),
            #('experiment_cVDPV2_NGA_100km_baseline_RI_SIA02N', 9),
            ]

# *****************************************************************************


def make_fig():

    dy_init = 1
    dy_end = 0
    c_thresh = 1

    tpath = os.path.join('..', 'Assets', 'data','routine_dat.json')
    with open(tpath) as fid01:
        ri_rate = json.load(fid01)

    tpath = os.path.join('..', 'Assets', 'data','shapes_adm00.json')
    with open(tpath) as fid01:
        adm00_shp = json.load(fid01)

    tpath = os.path.join('..', 'Assets', 'data','shapes_adm02.json')
    with open(tpath) as fid01:
        adm02_shp = json.load(fid01)

    tpath = os.path.join('..', 'Assets', 'data','epi_dat_mo.json')
    with open(tpath) as fid01:
        epi_dat_mo = json.load(fid01)

    for dir_tup in DIRNAMES:
 
        dirname = dir_tup[0]
        fig_clr = 'C{:d}'.format(dir_tup[1])

        # Sim outputs
        tpath = os.path.join('..', dirname)

        with open(os.path.join(tpath, 'data_brick.json')) as fid01:
            data_brick = json.load(fid01)

        with open(os.path.join(tpath, P_FILE)) as fid01:
            param_dict = json.load(fid01)

        START_YEAR = int(param_dict[EXP_C]['start_year'])
        NUM_YEARS = int(param_dict[EXP_C]['run_years'])
        N_SIMS = param_dict[NUM_SIMS]
        DT = int(t_step_days)

        NODE_DICT = data_brick['node_names']

        year_init = START_YEAR + dy_init
        run_years = NUM_YEARS - dy_init - dy_end

        MBIN = np.cumsum([0] + NUM_YEARS*MO_DAYS)
        t_vec = START_YEAR + np.arange(0, NUM_YEARS, 1/12) + 1/24
        tbool = (t_vec >= year_init) & (t_vec < (year_init+run_years))
        inf_data = np.zeros((N_SIMS, len(NODE_DICT), 12*NUM_YEARS))
        cal_data = np.zeros((N_SIMS, 2))

        afp_rate = np.zeros(len(NODE_DICT), dtype=float)
        for n_name in NODE_DICT:
            afp_rate[NODE_DICT[n_name]] = (1.0-ri_rate[n_name]*0.50)/850.0

        tvec_ref = np.array(epi_dat_mo[targ_adm00[0]]['times'])
        tbool_ref = (tvec_ref >= year_init) & (tvec_ref < (year_init+run_years))
        tvec_ref = tvec_ref[tbool_ref]
        ref_dat_mo = dict()
        ref_dat_mo_tot = None
        for cname in targ_adm00:
            ref_dat_mo[cname] = np.array(epi_dat_mo[cname]['cases'])[tbool_ref]
            if ref_dat_mo_tot is not None:
                ref_dat_mo_tot += ref_dat_mo[cname]
            else:
                ref_dat_mo_tot = ref_dat_mo[cname]

        ref_adm_idx = dict()
        for cname in targ_adm00:
            c_idx = np.array([NODE_DICT[adm02] for adm02 in NODE_DICT if adm02.startswith(cname+':')])
            n_idx = np.zeros(len(NODE_DICT), dtype=bool)
            n_idx[c_idx] = 1
            ref_adm_idx[cname] = n_idx

        for sim_idx_str in data_brick:
            if (not sim_idx_str.isdecimal()):
                continue
            sim_idx = int(sim_idx_str)
            inf_mat = np.array(data_brick[sim_idx_str]['infmat'])[:, :-1]
            inf_days = np.zeros((len(NODE_DICT), 365*NUM_YEARS))

            for k1 in range(DT):
                inf_days[:, k1::DT] = inf_mat/t_step_days

            for k2 in range(len(MBIN)-1):
                inf_data[sim_idx, :, k2] = np.sum(inf_days[:, MBIN[k2]:MBIN[k2+1]], axis=1)

            for cname in targ_adm00:
                tot_inf = np.sum(inf_data[sim_idx, ref_adm_idx[cname], :], axis=0).tolist()
                cal_tuple = norpois_opt(ref_dat_mo[cname], tot_inf)
                cal_data[sim_idx, 0] += cal_tuple[0]
                cal_data[sim_idx, 1] = 1/cal_tuple[1]

        totinf = np.sum(inf_data, axis=1)
        cuminf = np.cumsum(totinf, axis=1)
        gidx = (cuminf[:, -1] >= init_ob_thresh)
        #gidx = gidx & (cuminf[:, -1] > 900e3) #& (cuminf[:, -1] < 180e3)
        #gidx = gidx & (cuminf[:, -1] > 150e3) & (cuminf[:, -1] < 280e3)
        #gidx = gidx & (cuminf[:, -73] < 70e3)
        #gidx = gidx & (np.max(totinf, axis=1) < 6e3)
        #gidx = gidx & (totinf[:, -1] > 0) #& (cuminf[:, -1] < 180e3) #& (cuminf[:, -1] > 120e3)
        #gidx = gidx & (np.array(list(range(N_SIMS))) == 14) #& (np.array(list(range(N_SIMS))) < 900) #104
        #gidx = gidx & (np.array(list(range(N_SIMS))) > 200) #& (np.array(list(range(N_SIMS))) <= 225)

        print(np.sum(gidx))
        #print(np.argwhere(gidx), cuminf[gidx, -1])

        if (False):
            cal_list = np.argsort(cal_data[:,0])[::-1]
            idx_ge = 0
            for k1 in range(cal_list.shape[0]):
                if (gidx[cal_list[k1]]):
                    print(cal_data[cal_list[k1],:], cal_list[k1], cuminf[cal_list[k1], -1])
                    idx_ge += 1
                if (idx_ge > 10):
                    break

        #print(np.sum(gidx))
        #print(np.argwhere(gidx))
        #gidx = (np.array(list(range(N_SIMS))) == 652)

        if (False):
            for n7 in range(gidx.shape[0]):
                idx_str = '{:05d}'.format(n7)
                if (not gidx[n7] and idx_str in data_brick):
                    data_brick.pop('{:05d}'.format(n7))
            with open(os.path.join(tpath, 'data_brick_2.json'), 'w') as fid01:
                json.dump(data_brick, fid01)
            1/0

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

        #obp_lab = 'Fraction: {:5.3f}'.format(np.sum(gidx)/N_SIMS)
        #axs01.text(0.05, 0.9, obp_lab, fontsize=14, transform = axs01.transAxes)

        yval1 = totinf[gidx]/1000
        yval2 = np.mean(yval1, axis=0)
        for k3 in range(yval1.shape[0]):
            #axs01.plot(t_vec[tbool], yval1[k3, tbool], '.', c=fig_clr, alpha=0.1)
            axs01.plot(t_vec[tbool], yval1[k3, tbool])
        #axs01.plot(t_vec[tbool], yval2[tbool], c='k', lw=3)


        axs01.set_xlim(year_init, year_init+run_years+0.002)
        ticloc01 = list(range(year_init, year_init+run_years+1))
        ticlab01 = [str(val) for val in ticloc01]
        axs01.set_xticks(ticks=ticloc01)
        axs01.set_xticklabels(ticlab01)

        axs01.set_ylabel('Simulated Incidence (thousands)', fontsize=18)
        axs01.set_yscale('symlog', linthresh=100)
        axs01.set_ylim(0, 400)
        ticloc02 = list(range(0,241,40))+list(range(320,801,80))
        ticlab02 = [str(val) for val in ticloc02]
        ticlab02 = ticlab02[:9]+4*['']+[ticlab02[-1]]
        axs01.set_yticks(ticks=ticloc02)
        #axs01.set_yticklabels(ticlab02)
        axs01.set_yticklabels(len(ticloc02)*[''])

        axs02 = axs01.twinx()
        axs02.bar(tvec_ref, ref_dat_mo_tot, width=1/12,
                  alpha=0.2, facecolor='C3', edgecolor=None)

        #totafp1 = np.sum(inf_data*afp_rate[np.newaxis, :, np.newaxis], axis=1)
        #totafp2 = np.mean(totafp1[gidx], axis=0)
        #totafp2a = np.sum(totafp1[gidx,(8*73):(13*73)], axis=1)
        #print(dirname)
        #print(np.sum(totafp2[:550]))
        #print(np.sum(ref_dat_mo_tot))
        #print(np.mean(totafp2a))
        #print(np.quantile(totafp2a,[0.05, 0.95]))
        #print()

        #totafp3 = totafp2*73/12 # Time step conversion for visualization
        #axs02.plot(t_vec[tbool], totafp3[tbool], c='C3', lw=2, ls='--')

        axs02.set_yscale('symlog', linthresh=100)
        axs02.set_ylim(0, 400)
        ticloc02 = list(range(0,121,20))+list(range(160,401,40))
        ticlab02 = [str(val) for val in ticloc02]
        ticlab02 = ticlab02[:9]+4*['']+[ticlab02[-1]]
        axs02.set_yticks(ticks=ticloc02)
        axs02.set_yticklabels(ticlab02)
        axs02.tick_params(axis='y', which='major', labelsize=14)
        axs02.set_ylabel('Observed AFP Cases', fontsize=18, c='C3')

        (xmin, xmax, ymin, ymax) = (180, -180, 90, -90)
        for cname in targ_adm00:
            pts = np.array(adm00_shp[cname]['points'])
            xmin = np.floor(min(xmin, np.min(pts[:,0])))
            xmax = np.ceil(max(xmax, np.max(pts[:,0])))
            ymin = np.floor(min(ymin, np.min(pts[:,1])))
            ymax = np.ceil(max(ymax, np.max(pts[:,1])))
        xydelt = max(xmax-xmin, ymax-ymin)

        adm02_mat01 = (inf_data[gidx,:,:]>c_thresh)
        adm02_mat02 = (inf_data[gidx,:,:]>1)
        for k1 in range(run_years):
            axs01 = axlist[k1+1]
            axs01.axis('off')
            axs01.set_xlim(xmin, xmin+xydelt)
            axs01.set_ylim(ymin, xmin+xydelt)
            for cname in targ_adm00:
                adm00_prt = adm00_shp[cname]['parts']
                adm00_pts = adm00_shp[cname]['points']
                shape_patch(axs01, adm00_pts, adm00_prt, clr=3*[0.9])
                shape_line(axs01, adm00_pts, adm00_prt, wid=0.2)
            yidx = (t_vec>=ticloc01[k1]) & (t_vec<ticloc01[k1+1])
            yrdat01 = np.max(adm02_mat01[:, :, yidx], axis=2)
            yrdat01 = np.mean(yrdat01, axis=0)
            yrdat02 = np.max(adm02_mat02[:, :, yidx], axis=2)
            yrdat02 = np.mean(yrdat02, axis=0)
            for adm02_name in NODE_DICT:
                k2 = NODE_DICT[adm02_name]
                if (yrdat01[k2] > 0):
                    adm02_prt = adm02_shp[adm02_name]['parts']
                    adm02_pts = adm02_shp[adm02_name]['points']
                    clr_val = [1.0, 1.0-yrdat01[k2], 1.0-yrdat01[k2]]
                    shape_patch(axs01, adm02_pts, adm02_prt, clr=clr_val)
                elif (yrdat02[k2] > 0):
                    adm02_prt = adm02_shp[adm02_name]['parts']
                    adm02_pts = adm02_shp[adm02_name]['points']
                    clr_val = [1.0, 1.0, 1.0]
                    shape_patch(axs01, adm02_pts, adm02_prt, clr=clr_val)

        plt.tight_layout()
        plt.savefig('fig_extent_{:s}_01_v5.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
