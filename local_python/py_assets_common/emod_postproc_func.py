# *****************************************************************************
#
# *****************************************************************************

import json
import os
import struct

import numpy as np

from emod_constants import AGE_KEY_LIST, YR_DAYS, POP_PYR, CBR_VEC, \
                           NODE_IDS_STR, NODE_POP_STR, INF_FRAC, RST_FILE, \
                           RST_TIME, RST_CONT_INF, RST_CONT_TOT, R0_VEC, \
                           R0_TIME, CAMP_COST

# *****************************************************************************


def post_proc_poppyr(output_path, parsed_out):

    # Sample population pyramid every year
    with open(os.path.join(output_path, 'DemographicsSummary.json')) as fid01:
        demog_output = json.load(fid01)

    ds_start = demog_output['Header']['Start_Time']
    ds_nstep = demog_output['Header']['Timesteps']
    ds_tsize = demog_output['Header']['Simulation_Timestep']
    time_vec = np.arange(ds_nstep)*ds_tsize + ds_start
    nyr_bool = (np.diff(time_vec//YR_DAYS) > 0.0)
    run_yrs = ds_nstep*ds_tsize/YR_DAYS

    pyr_dat = np.zeros((int(run_yrs)+1, len(AGE_KEY_LIST)))

    for k1 in range(len(AGE_KEY_LIST)):
        age_key_str = 'Population Age {:s}'.format(AGE_KEY_LIST[k1])
        age_vec_dat = np.array(demog_output['Channels'][age_key_str]['Data'])
        pyr_dat[0, k1] = age_vec_dat[0]
        age_subset = age_vec_dat[:-1][nyr_bool]
        if (age_subset.shape[0] < int(run_yrs)):
            age_subset = np.append(age_subset, age_vec_dat[-1])
        pyr_dat[1:, k1] = age_subset

    parsed_out[POP_PYR] = pyr_dat.tolist()

    return None

# *****************************************************************************


def post_proc_cbr(output_path, parsed_out):

    # Retain annualized count of births
    with open(os.path.join(output_path, 'InsetChart.json')) as fid01:
        inset_chart = json.load(fid01)

    ic_start = inset_chart['Header']['Start_Time']
    ic_nstep = inset_chart['Header']['Timesteps']
    ic_tsize = inset_chart['Header']['Simulation_Timestep']
    cumb_vec = np.array(inset_chart['Channels']['Births']['Data'])

    time_vec = np.arange(ic_nstep)*ic_tsize + ic_start
    nyr_bool = (np.diff(time_vec//365.0) > 0.0)
    run_years = (ic_nstep*ic_tsize)//365.0
    b_vec = cumb_vec[:-1][nyr_bool]
    if (b_vec.shape[0] < run_years):
        b_vec = np.append(b_vec, cumb_vec[-1])
    b_vec[1:] = np.diff(b_vec)

    parsed_out[CBR_VEC] = b_vec.tolist()

    return None

# *****************************************************************************


def post_proc_nodepop(output_path, parsed_out):

    # Population data from Spatial_Output_Channels = ["Population"]
    fname = 'SpatialReport_Population.bin'
    with open(os.path.join(output_path, fname), mode='rb') as fid01:
        sr_data = fid01.read()

    # Struct output is tuple even when single value
    num_nodes = struct.unpack("i", sr_data[0:4])[0]
    num_times = struct.unpack("i", sr_data[4:8])[0]
    node_ids = struct.unpack("i"*num_nodes,
                             sr_data[8:(8+4*num_nodes)])
    sim_data = struct.unpack("f"*num_nodes*num_times,
                             sr_data[(8+4*num_nodes):])

    # Construct numpy data structures
    node_id_vec = np.array([val for val in node_ids])
    pop_mat = np.reshape(np.array([val for val in sim_data]), (num_times, -1))

    parsed_out[NODE_IDS_STR] = node_id_vec.tolist()
    parsed_out[NODE_POP_STR] = pop_mat.tolist()

    return None

# *****************************************************************************


def post_proc_prev(output_path, parsed_out):

    # Retain timeseries of infected fraction
    with open(os.path.join(output_path, 'InsetChart.json')) as fid01:
        inset_chart = json.load(fid01)

    inf_frac_vec = np.array(inset_chart['Channels']['Infected']['Data'])

    parsed_out[INF_FRAC] = inf_frac_vec.tolist()

    return None

# *****************************************************************************


def post_proc_cost(output_path, parsed_out):

    # Retain timeseries of campaign cost
    with open(os.path.join(output_path, 'InsetChart.json')) as fid01:
        inset_chart = json.load(fid01)

    inf_frac_vec = np.array(inset_chart['Channels']['Campaign Cost']['Data'])

    parsed_out[CAMP_COST] = inf_frac_vec.tolist()

    return None

# *****************************************************************************


def post_proc_R0(output_path, parsed_out):

    # Retain timeseries of infected fraction
    with open(os.path.join(output_path, RST_FILE)) as fid01:
        rst_dat = np.loadtxt(fid01, delimiter=',', skiprows=1)

    rst_time_vec = rst_dat[:, RST_TIME]
    rst_cont_vec = rst_dat[:, RST_CONT_TOT]
    rst_infs_vec = rst_dat[:, RST_CONT_INF]
    np_eps = np.finfo(float).eps

    tot_contagion = rst_cont_vec
    tot_infection = rst_infs_vec

    est_r0_vec = tot_contagion/(tot_infection+np_eps)

    parsed_out[R0_TIME] = rst_time_vec.tolist()
    parsed_out[R0_VEC] = est_r0_vec.tolist()

    return None

# *****************************************************************************
