# *****************************************************************************

import json
import os
import sys

import numpy as np

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))
from py_assets_common.emod_local_proc import crs_proc
from py_assets_common.emod_constants import EXP_C, EXP_V, CBR_VEC, \
                                            NUM_SIMS, P_FILE, POP_PYR
from global_data import run_years, start_year
from ref_dat import pop_2019

# *****************************************************************************


def make_dat():

    with open('data_brick.json') as fid01:
        data_brick = json.load(fid01)

    with open(P_FILE) as fid01:
        param_dict = json.load(fid01)

    nsims = int(param_dict[NUM_SIMS])
    ss_demog = param_dict[EXP_C]['steady_state_demog']
    demog_set = param_dict[EXP_C]['demog_set']
    adm01 = [param_dict[EXP_C]['adm01']]

    prov_list = sorted(list(set(adm01)))

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

    for prov in prov_list:
        gidx = (np.array(adm01) == prov) & fidx

        # Average population
        pyr_mat_avg = np.mean(pyr_mat[gidx, :, :], axis=0)
        pop_tot = np.sum(pyr_mat_avg, axis=1)
        pop_tot = np.diff(pop_tot)/2.0 + pop_tot[:-1]

        # CRS calculations
        XDAT = np.arange(start_year, start_year+run_years) + 0.5
        fname = 'fert_dat_{:s}.csv'.format(demog_set)
        fnabs = os.path.abspath(os.path.join('..', 'Assets', 'data', fname))
        (frt_brth, crs_prob_vec) = crs_proc(fnabs, XDAT, pyr_mat_avg, ss_demog)

        # Normalize timeseries required for CRS calculation
        brth_vec = np.mean(birth_mat[gidx, :], axis=0)
        norm_crs_timevec = brth_vec/frt_brth
        print(prov, np.mean(norm_crs_timevec))

        # Population scaling
        pop_scale = np.sum([pop_2019[val] for val in pop_2019])/pop_tot[19]

        # Calculation CRS from infections and probabilities
        crs_mat = np.sum(inf_mat[gidx, :, :]*crs_prob_vec[:,0], axis=2)
        crs_out = crs_mat*pop_scale

        # CRS range
        crs_out = np.cumsum(crs_out[:, -35:-5], axis=1)
        crs_sort = np.sort(crs_out, axis=1)
        tidx = int(crs_sort.shape[0]/20)
        print(prov, int(np.mean(crs_sort[:,-1])))

        # CRS rates
        ydat = np.mean(crs_mat, axis=0)/brth_vec*1e5
        print(prov, np.round(np.mean(ydat[-15:-5])))

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_dat()

# *****************************************************************************