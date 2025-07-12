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
                                            CAMP_COST, MO_DAYS

from py_assets_common.emod_analysis import norpois_opt

from global_data import base_year, init_ob_thresh, targ_adm00, t_step_days

# *****************************************************************************

DIRNAMES = [
            ('experiment_cVDPV2_100km_base', 0),
            ('experiment_cVDPV2_100km_base_ri2027', 0),
            ('experiment_cVDPV2_100km_base_obr2026', 0),
            ('experiment_cVDPV2_100km_base_obr2026_ri2027', 0),

            ('experiment_cVDPV2_100km_base_sia01-NGA', 1),
            ('experiment_cVDPV2_100km_base_sia01-NGA_ri2027', 0),
            ('experiment_cVDPV2_100km_base_sia01u2-NGA', 0),
            ('experiment_cVDPV2_100km_base_sia01u2-NGA_ri2027', 0),
            ('experiment_cVDPV2_100km_base_sia02-NGA', 0),
            ('experiment_cVDPV2_100km_base_sia02-NGA_ri2027', 0),
            ('experiment_cVDPV2_100km_base_sia02u2-NGA', 0),
            ('experiment_cVDPV2_100km_base_sia02u2-NGA_ri2027', 0),

            ('experiment_cVDPV2_100km_base_sia01-NGAN', 0),
            ('experiment_cVDPV2_100km_base_sia01-NGAN_ri2027', 0),
            ('experiment_cVDPV2_100km_base_sia01u2-NGAN', 0),
            ('experiment_cVDPV2_100km_base_sia01u2-NGAN_ri2027', 0),
            ('experiment_cVDPV2_100km_base_sia02-NGAN', 0),
            ('experiment_cVDPV2_100km_base_sia02-NGAN_ri2027', 0),
            ('experiment_cVDPV2_100km_base_sia02u2-NGAN', 0),
            ('experiment_cVDPV2_100km_base_sia02u2-NGAN_ri2027', 0),

            ('experiment_cVDPV2_100km_base_sia01-LKCHAD', 0),
            ('experiment_cVDPV2_100km_base_sia01-LKCHAD_ri2027', 0),
            ('experiment_cVDPV2_100km_base_sia01u2-LKCHAD', 0),
            ('experiment_cVDPV2_100km_base_sia01u2-LKCHAD_ri2027', 0),
            ('experiment_cVDPV2_100km_base_sia02-LKCHAD', 0),
            ('experiment_cVDPV2_100km_base_sia02-LKCHAD_ri2027', 0),
            ('experiment_cVDPV2_100km_base_sia02u2-LKCHAD', 0),
            ('experiment_cVDPV2_100km_base_sia02u2-LKCHAD_ri2027', 0),

            #('experiment_cVDPV2_100km_base_sia1p1-LKCHAD_obr2026', 0),
            #('experiment_cVDPV2_100km_base_sia1p1-NGA_obr2026', 0),
            #('experiment_cVDPV2_100km_base_sia1p1-NGAN_obr2026', 0),

            ('experiment_cVDPV2_100km_base_sia2p1-LKCHAD', 0),
            #('experiment_cVDPV2_100km_base_sia2p1-LKCHAD_obr2026', 0),
            ('experiment_cVDPV2_100km_base_sia2p1-NGA', 0),
            #('experiment_cVDPV2_100km_base_sia2p1-NGA_obr2026', 0),
            ('experiment_cVDPV2_100km_base_sia2p1-NGAN', 0),
            #('experiment_cVDPV2_100km_base_sia2p1-NGAN_obr2026', 0),

            ('experiment_cVDPV2_100km_base_sia3p1-LKCHAD', 0),
            #('experiment_cVDPV2_100km_base_sia3p1-LKCHAD_obr2026', 0),
            ('experiment_cVDPV2_100km_base_sia3p1-NGA', 0),
            #('experiment_cVDPV2_100km_base_sia3p1-NGA_obr2026', 0),
            ('experiment_cVDPV2_100km_base_sia3p1-NGAN', 0),
            #('experiment_cVDPV2_100km_base_sia3p1-NGAN_obr2026', 0),

            ('experiment_cVDPV2_100km_base_sia4p1-LKCHAD', 0),
            #('experiment_cVDPV2_100km_base_sia4p1-LKCHAD_obr2026', 0),
            ('experiment_cVDPV2_100km_base_sia4p1-NGA', 0),
            #('experiment_cVDPV2_100km_base_sia4p1-NGA_obr2026', 0),
            ('experiment_cVDPV2_100km_base_sia4p1-NGAN', 0),
            #('experiment_cVDPV2_100km_base_sia4p1-NGAN_obr2026', 0),
            ]

# *****************************************************************************


def make_fig():

    dy_init = 3
    dy_end = 0

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

        NODE_DICT = data_brick['node_names']
        START_YEAR = int(param_dict[EXP_C]['start_year'])
        NUM_YEARS = int(param_dict[EXP_C]['run_years'])
        N_SIMS = param_dict[NUM_SIMS]
        DT = int(t_step_days)
        MBIN = np.cumsum([0] + NUM_YEARS*MO_DAYS)

        year_init = START_YEAR + dy_init
        run_years = NUM_YEARS - dy_init - dy_end

        t_vec = np.array(data_brick['t_vec'])/365 + base_year
        tbool = (t_vec >= year_init) & (t_vec < (year_init+run_years))
        inf_data = np.zeros((N_SIMS, len(NODE_DICT), 12*NUM_YEARS))
        cal_data = np.zeros((N_SIMS, 2))
        cal_data[:, 1] = 1
        sia_data = np.zeros((N_SIMS, t_vec.shape[0]))

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
            siavec = np.array(data_brick[sim_idx_str][CAMP_COST])
            sia_data[sim_idx, :] = siavec
            for k1 in range(DT):
                inf_days[:, k1::DT] = inf_mat/t_step_days

            for k2 in range(len(MBIN)-1):
                inf_data[sim_idx, :, k2] = np.sum(inf_days[:, MBIN[k2]:MBIN[k2+1]], axis=1)

            tot_inf = np.sum(inf_data[sim_idx, :, :], axis=0).tolist()
            cal_tuple = norpois_opt(ref_dat_mo_tot, tot_inf)
            cal_data[sim_idx, 0] = cal_tuple[0]
            cal_data[sim_idx, 1] = 984.664 #1/cal_tuple[1]

        tot_inf = np.sum(inf_data, axis=1)
        cum_inf = np.cumsum(tot_inf, axis=1)
        gidx0 = (cum_inf[:, -1] >= init_ob_thresh)

        #gidx = gidx & (cum_inf[:, -1] > 900e3) #& (cum_inf[:, -1] < 180e3)
        #gidx = gidx & (cum_inf[:, -1] > 150e3)
        #gidx = gidx & (np.sum(tot_inf[:, -96:-90], axis=1) > 0)
        #gidx = gidx & (np.array(list(range(n_sims))) == 55) #& (np.array(list(range(n_sims))) < 900) #104

        #print(np.sum(gidx))
        #print(np.argwhere(gidx))
        #print(cum_inf[gidx, -1])

        gidx = gidx0 & (np.sum(tot_inf[:, -12:], axis=1) > 0)

        #dcases = (cum_inf[:, -1] - cum_inf[:, -60])/1200
        #dcidx = np.argsort(dcases)[-200:]
        ddoses = (sia_data[gidx, -1] - sia_data[gidx, -(73*5)])/1e6
        #mdoses = ddoses[dcidx]
        mean_val = np.mean(ddoses)
        quant_val = np.quantile(ddoses, [0.05, 0.95])

        str_out = ''
        str_out = str_out + str(np.sum(gidx)) + ','
        str_out = str_out + dirname + ','
        str_out = str_out + str(int(quant_val[0])) + ','
        str_out = str_out + str(int(mean_val)) + ','
        str_out = str_out + str(int(quant_val[1]))

        print(str_out)

        continue

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
        #axs01.set_ylim(0, 120)

        plt.tight_layout()
        plt.savefig('fig_doses_{:s}_02.png'.format(dirname))
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
