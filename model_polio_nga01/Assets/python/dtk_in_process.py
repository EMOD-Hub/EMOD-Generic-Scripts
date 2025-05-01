# *****************************************************************************
#
# *****************************************************************************

import os

import global_data as gdata

import numpy as np

from emod_constants import RST_FILE, RST_TIME, RST_NODE, RST_CLADE, \
                           RST_GENOME, RST_TOT_INF, RST_NEW_INF


# *****************************************************************************


def application(timestep):

    TIME_VAL = float(timestep)
    T_DELTA = gdata.inproc_dt
    T_STEP = gdata.t_step_days
    T_START = TIME_VAL-T_DELTA

    # Variables for this simulation
    USE_OBR = gdata.var_params['use_obr']

    # Evaluate outbreak status every dt days
    if ( not USE_OBR or (TIME_VAL%365)%T_DELTA ):
        return None

    # Load strain data
    infdat = np.loadtxt(os.path.join('output', RST_FILE),
                         delimiter=',', skiprows=1, ndmin=2)

    if (infdat.shape[0] == 0):
        return None

    ADM01_DICT = gdata.adm01_idlist
    ADM01_LIST = sorted(list(ADM01_DICT.keys()))
    ADM01_IDX = {val[0]: val[1] for val in
                 zip(ADM01_LIST, range(len(ADM01_LIST)))}

    # Construct array of cVDPV cases
    cVDPV_genome = gdata.boxes_sabin2 + gdata.boxes_nopv2
    ntstep = int(T_DELTA/T_STEP)+1
    dbrick0 = np.zeros((len(ADM01_LIST), ntstep), dtype=int)

    for adm01 in ADM01_LIST:
        odex = ADM01_IDX[adm01]
        gidx = (infdat[:, RST_TIME] >= T_START)
        gidx = (infdat[:, RST_CLADE] == 0) & gidx
        gidx = (infdat[:, RST_GENOME] == cVDPV_genome) & gidx
        gidx = np.isin(infdat[:, RST_NODE], ADM01_DICT[adm01]) & gidx
        subdat = infdat[gidx, :]
        tvec = (subdat[:, RST_TIME]-T_START)/T_STEP
        for k1 in range(subdat.shape[0]):
            tdex = int(tvec[k1])
            dbrick0[odex, tdex] += subdat[k1, RST_NEW_INF]

    print(np.sum(dbrick0, axis=1))
    dbrick1 = np.random.binomial(dbrick0, 0.005)
    print(np.sum(dbrick1, axis=1))


    return None

# *****************************************************************************
