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

from py_assets_common.emod_analysis import norpois_opt
from py_assets_common.emod_local_proc import shape_patch, shape_line

from global_data import base_year, init_ob_thresh, targ_adm00, t_step_days

# *****************************************************************************

DIRNAMES = [
            #('experiment_cVDPV2_100km_base', 0),
            #('experiment_cVDPV2_100km_base_ri2027', 0),
            #('experiment_cVDPV2_100km_base_obr2026', 0),
            #('experiment_cVDPV2_100km_base_obr2026_ri2027', 0),

            ('experiment_cVDPV2_100km_base_obr2026y1_ri2027d2', 4),
            ('experiment_cVDPV2_100km_base_sia01p1-NGAN_ri2027d2', 5),

            ('experiment_cVDPV2_100km_base_sia01-NGA', 1),


            #('experiment_cVDPV2_100km_base_sia01-NGA_ri2027', 0),
            #('experiment_cVDPV2_100km_base_sia01u2-NGA', 0),
            #('experiment_cVDPV2_100km_base_sia01u2-NGA_ri2027', 0),
            #('experiment_cVDPV2_100km_base_sia02-NGA', 0),
            #('experiment_cVDPV2_100km_base_sia02-NGA_ri2027', 0),
            #('experiment_cVDPV2_100km_base_sia02u2-NGA', 0),
            #('experiment_cVDPV2_100km_base_sia02u2-NGA_ri2027', 0),

            #('experiment_cVDPV2_100km_base_sia01-NGAN', 0),
            #('experiment_cVDPV2_100km_base_sia01-NGAN_ri2027', 0),
            #('experiment_cVDPV2_100km_base_sia01u2-NGAN', 0),
            #('experiment_cVDPV2_100km_base_sia01u2-NGAN_ri2027', 0),
            #('experiment_cVDPV2_100km_base_sia02-NGAN', 0),
            #('experiment_cVDPV2_100km_base_sia02-NGAN_ri2027', 0),
            #('experiment_cVDPV2_100km_base_sia02u2-NGAN', 0),
            #('experiment_cVDPV2_100km_base_sia02u2-NGAN_ri2027', 0),

            #('experiment_cVDPV2_100km_base_sia01-LKCHAD', 0),
            #('experiment_cVDPV2_100km_base_sia01-LKCHAD_ri2027', 0),
            #('experiment_cVDPV2_100km_base_sia01u2-LKCHAD', 0),
            #('experiment_cVDPV2_100km_base_sia01u2-LKCHAD_ri2027', 0),
            #('experiment_cVDPV2_100km_base_sia02-LKCHAD', 0),
            #('experiment_cVDPV2_100km_base_sia02-LKCHAD_ri2027', 0),
            #('experiment_cVDPV2_100km_base_sia02u2-LKCHAD', 0),
            #('experiment_cVDPV2_100km_base_sia02u2-LKCHAD_ri2027', 0),
            ]

# *****************************************************************************


def make_fig():

    dy_init = 3
    dy_end = 0
    c_thresh = 1

    #tpath = os.path.join('..', 'Assets', 'data','routine_dat.json')
    #with open(tpath) as fid01:
    #    ri_rate = json.load(fid01)

    tpath = os.path.join('..', 'Assets', 'data','shapes_adm00.json')
    with open(tpath) as fid01:
        adm00_shp = json.load(fid01)

    tpath = os.path.join('..', 'Assets', 'data','shapes_adm02.json')
    with open(tpath) as fid01:
        adm02_shp = json.load(fid01)

    tpath = os.path.join('..', 'Assets', 'data','epi_dat_mo.json')
    with open(tpath) as fid01:
        epi_dat_mo = json.load(fid01)

    loopone = True
    for dir_tup in DIRNAMES:
 
        dirname = dir_tup[0]
        fig_clr = 'C{:d}'.format(dir_tup[1])

        # Sim outputs
        tpath = os.path.join('..', dirname)

        with open(os.path.join(tpath, 'data_brick.json')) as fid01:
            data_brick = json.load(fid01)

        with open(os.path.join(tpath, P_FILE)) as fid01:
            param_dict = json.load(fid01)

        NODE_DICT = data_brick['node_names']
        START_YEAR = int(param_dict[EXP_C]['start_year'])
        NUM_YEARS = int(param_dict[EXP_C]['run_years'])
        N_SIMS = param_dict[NUM_SIMS]
        DT = int(t_step_days)
        MBIN = np.cumsum([0] + NUM_YEARS*MO_DAYS)

        year_init = START_YEAR + dy_init
        run_years = NUM_YEARS - dy_init - dy_end

        t_vec = START_YEAR + np.arange(0, NUM_YEARS, 1/12) + 1/24
        tbool = (t_vec >= year_init) & (t_vec < (year_init+run_years))
        inf_data = np.zeros((N_SIMS, len(NODE_DICT), 12*NUM_YEARS))
        cal_data = np.zeros((N_SIMS, 2))
        cal_data[:, 1] = 1

        #afp_rate = np.zeros(len(NODE_DICT), dtype=float)
        #for n_name in NODE_DICT:
        #    afp_rate[NODE_DICT[n_name]] = (1.0-ri_rate[n_name]*0.50)/850.0

        tvec_ref = np.array(epi_dat_mo[targ_adm00[0]]['times'])
        tbool_ref = (tvec_ref >= year_init) & (tvec_ref < (year_init+run_years))
        tvec_ref = tvec_ref[tbool_ref]
        ref_dat_mo_tot = None
        for cname in targ_adm00:
            ref_dat_mo = np.array(epi_dat_mo[cname]['cases'])[tbool_ref]
            if ref_dat_mo_tot is not None:
                ref_dat_mo_tot += ref_dat_mo
            else:
                ref_dat_mo_tot = ref_dat_mo

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

            tot_inf = np.sum(inf_data[sim_idx, :, :], axis=0).tolist()
            cal_tuple = norpois_opt(ref_dat_mo_tot, tot_inf)
            cal_data[sim_idx, 0] = cal_tuple[0]
            cal_data[sim_idx, 1] = 984.664 #1/cal_tuple[1]
            #if(np.sum(tot_inf[-24:-12]) < 200e3):
            #   cal_data[sim_idx, 0] += 1000 
            #if(np.sum(tot_inf[-12:]) > 100e3):
            #   cal_data[sim_idx, 0] += 1000 

        tot_inf = np.sum(inf_data, axis=1)
        cum_inf = np.cumsum(tot_inf, axis=1)
        gidx = (cum_inf[:, -1] >= init_ob_thresh)
        #gidx = gidx & (cum_inf[:, -1] > 900e3) #& (cum_inf[:, -1] < 180e3)
        #gidx = gidx & (cum_inf[:, -1] > 150e3) & (cum_inf[:, -1] < 280e3)
        #gidx = gidx & (np.sum(tot_inf[:, -12:], axis=1) > 0)
        #gidx = gidx & (np.array(list(range(N_SIMS))) == 14) #& (np.array(list(range(N_SIMS))) < 900)

        #print(np.sum(gidx))
        #print(np.argwhere(gidx), cum_inf[gidx, -1])

        if (False):
            cal_list = np.argsort(cal_data[:,0])[::-1]
            idx_ge = 0
            for k1 in range(cal_list.shape[0]):
                if (gidx[cal_list[k1]]):
                    print(cal_data[cal_list[k1],:], cal_list[k1], cum_inf[cal_list[k1], -1])
                    idx_ge += 1
                if (idx_ge > 10):
                    break

        #print(np.sum(gidx))
        #print(np.argwhere(gidx))
        #gidx = (np.array(list(range(N_SIMS))) == 249)

        dcases = (cum_inf[:, -1] - cum_inf[:, -60])/1200
        dcases = np.sort(dcases)[-200:]
        mean_val = np.mean(dcases)
        quant_val = np.quantile(dcases, [0.05, 0.95])

        str_out = ''
        str_out = str_out + str(np.sum(gidx)) + ','
        str_out = str_out + dirname + ','
        str_out = str_out + str(int(quant_val[0])) + ','
        str_out = str_out + str(int(mean_val)) + ','
        str_out = str_out + str(int(quant_val[1]))

        print(str_out)

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

        # Timeseries
        yval1 = tot_inf[gidx]/cal_data[gidx,1][..., np.newaxis]
        yval2 = np.mean(yval1, axis=0)
        for k3 in range(yval1.shape[0]):
            axs01.plot(t_vec[tbool], yval1[k3, tbool], '.', c=fig_clr, alpha=0.2)
            #axs01.plot(t_vec[tbool], yval1[k3, tbool])
        axs01.plot(t_vec[tbool], yval2[tbool], c='k', lw=3)

        axs01.set_xlim(year_init, year_init+run_years+0.002)
        ticloc_x = list(range(year_init, year_init+run_years+1))
        ticlab_x = [str(val) for val in ticloc_x]
        axs01.set_xticks(ticks=ticloc_x)
        axs01.set_xticklabels(ticlab_x)

        axs01.set_ylabel('Observed AFP Cases', fontsize=18)
        axs01.set_yscale('symlog', linthresh=100)
        axs01.set_ylim(0, 300)
        ticloc_y = list(range(0,101,20))+list(range(200,301,100))
        ticlab_y = [str(val) for val in ticloc_y]
        #ticlab_y[-2] = ''
        #ticlab_y[-4] = ''
        #ticlab_y[-6] = ''
        axs01.set_yticks(ticks=ticloc_y)
        axs01.set_yticklabels(ticlab_y)
        axs01.tick_params(axis='y', which='major', labelsize=14)

        axs01.bar(tvec_ref, ref_dat_mo_tot, width=1/12,
                  alpha=0.2, facecolor='C3', edgecolor=None)

        # Incidence maps
        (xmin, xmax, ymin, ymax) = (180, -180, 90, -90)
        for cname in targ_adm00:
            pts = np.array(adm00_shp[cname]['points'])
            xmin = np.floor(min(xmin, np.min(pts[:,0])))
            xmax = np.ceil(max(xmax, np.max(pts[:,0])))
            ymin = np.floor(min(ymin, np.min(pts[:,1])))
            ymax = np.ceil(max(ymax, np.max(pts[:,1])))
        xydelt = max(xmax-xmin, ymax-ymin)

        adm02_mat01 = (inf_data[gidx,:,:]>c_thresh)
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
            yidx = (t_vec>=ticloc_x[k1]) & (t_vec<ticloc_x[k1+1])
            yrdat01 = np.max(adm02_mat01[:, :, yidx], axis=2)
            yrdat01 = np.mean(yrdat01, axis=0)
            for adm02_name in NODE_DICT:
                k2 = NODE_DICT[adm02_name]
                if (yrdat01[k2] > 0):
                    adm02_prt = adm02_shp[adm02_name]['parts']
                    adm02_pts = adm02_shp[adm02_name]['points']
                    clr_val = [1.0, 1.0-yrdat01[k2], 1.0-yrdat01[k2]]
                    shape_patch(axs01, adm02_pts, adm02_prt, clr=clr_val)

        plt.tight_layout()
        plt.savefig('fig_extent_{:s}_01_v6.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
