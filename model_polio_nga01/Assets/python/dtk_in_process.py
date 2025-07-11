# *****************************************************************************
#
# *****************************************************************************

import json
import os

import global_data as gdata

import numpy as np

import emod_api.campaign as camp_module

from emod_camp_events import ce_OPV_SIA
from emod_constants import RST_FILE, RST_TIME, RST_NODE, RST_CLADE, \
                           RST_GENOME, RST_NEW_INF

# *****************************************************************************


def application(timestep):

    TIME_VAL = float(timestep)
    T_DELTA = gdata.inproc_dt
    T_LAST = TIME_VAL-T_DELTA
    T_STEP = gdata.t_step_days

    SIM_IDX = gdata.sim_index

    DELAY_DICT = gdata.inproc_dict_delay
    ADM01_DICT = gdata.adm01_idlist
    SIATIME_DICT = gdata.inproc_dict_sia_time
    SIA_BASE_TAKE = gdata.var_params['sia_base_vax_take']
    SIA_NOPV2_MOD = gdata.nopv2_sia_take_fac
    SIA_COVER = gdata.sia_cover_dict
    SIA_TAKE = SIA_BASE_TAKE*SIA_NOPV2_MOD

    OBR_START = gdata.var_params['obr_start']
    OBR_END = gdata.var_params['obr_end']
    DAY_OBR_START = 365.0*(OBR_START-gdata.base_year)
    DAY_OBR_END = 365.0*(OBR_END-gdata.base_year)

    # Evaluate outbreak status every dt days
    if (TIME_VAL < DAY_OBR_START or TIME_VAL > DAY_OBR_END):
        return None

    if ((TIME_VAL % 365) % T_DELTA):
        return None

    # Load delay paramters
    if (DELAY_DICT is None):
        fname = os.path.join('Assets', 'data', 'obr_lag_param.json')
        with open(fname) as fid01:
            DELAY_DICT = json.load(fid01)
        gdata.inproc_dict_delay = DELAY_DICT

    # Load time of most recent OBR SIA
    if (SIATIME_DICT is None):
        SIATIME_DICT = {adm01: 0.0 for adm01 in ADM01_DICT}
        gdata.inproc_dict_sia_time = SIATIME_DICT

    # Cycle random numbers
    np.random.rand(SIM_IDX)

    # Load strain data
    infdat = np.loadtxt(os.path.join('output', RST_FILE),
                        delimiter=',', skiprows=1, ndmin=2)

    if (infdat.shape[0] == 0):
        return None

    # Construct array of cVDPV cases
    ADM01_LIST = sorted(list(ADM01_DICT.keys()))
    cVDPV_genome = gdata.boxes_sabin2 + gdata.boxes_nopv2
    ntstep = int(T_DELTA/T_STEP)+1
    dbrick_inf = np.zeros((len(ADM01_LIST), ntstep), dtype=int)

    for k1 in range(len(ADM01_LIST)):
        adm01 = ADM01_LIST[k1]
        gidx = (infdat[:, RST_TIME] >= T_LAST)
        gidx = (infdat[:, RST_CLADE] == 0) & gidx
        gidx = (infdat[:, RST_GENOME] == cVDPV_genome) & gidx
        gidx = np.isin(infdat[:, RST_NODE], ADM01_DICT[adm01]) & gidx
        subdat = infdat[gidx, :]
        tvec = (subdat[:, RST_TIME]-T_LAST)/T_STEP
        for k2 in range(ntstep):
            dbrick_inf[k1, k2] = np.sum(subdat[(tvec == k2), RST_NEW_INF])

    dbrick_cases = np.random.binomial(dbrick_inf, 0.005)

    # New campaign file
    camp_module.reset()
    CAMP_FILE = 'campaign_{:05d}.json'.format(int(TIME_VAL))

    for k1 in range(len(ADM01_LIST)):
        if (np.sum(dbrick_cases[k1, :]) == 0):
            continue

        adm01 = ADM01_LIST[k1]
        adm00 = adm01.rsplit(':', 1)[0]
        shape = DELAY_DICT[adm00][0]
        scale = DELAY_DICT[adm00][1]
        min_day = SIATIME_DICT[adm01] + 180.0
        sia_day = None

        print(adm01, dbrick_cases[k1], TIME_VAL)

        for k2 in range(ntstep):
            cases = dbrick_cases[k1, k2]
            ctimes = np.random.gamma(shape, scale=scale, size=cases)
            ctimes = ctimes - T_STEP*(ntstep-k2) + TIME_VAL
            sia_targ = ctimes + 30.0
            sia_targ = sia_targ[sia_targ > min_day]
            if (sia_targ.shape[0]):
                new_day = np.min(sia_targ)
                if (sia_day):
                    sia_day = np.min([sia_day, new_day])
                else:
                    sia_day = new_day
        if (sia_day):
            SIATIME_DICT[adm01] = sia_day

        sia_tstep = SIATIME_DICT[adm01]
        if (sia_tstep > TIME_VAL and sia_tstep <= TIME_VAL+T_DELTA):
            sia_day = sia_tstep+T_STEP
            n_list = ADM01_DICT[adm01]
            n_dict = {nid: SIA_COVER[nid] for nid in n_list}
            camp_event = ce_OPV_SIA(n_dict, start_day=sia_day,
                                    take=SIA_TAKE, clade=1, genome=0,
                                    yrs_min=0.2, yrs_max=5.0)
            camp_module.add(camp_event)
            camp_event = ce_OPV_SIA(n_dict, start_day=sia_day+30.0,
                                    take=SIA_TAKE, clade=1, genome=0,
                                    yrs_min=0.2, yrs_max=5.0)
            camp_module.add(camp_event)
            print(adm01, TIME_VAL, sia_day)

    if (camp_module.campaign_dict["Events"]):
        camp_module.save(filename=CAMP_FILE)
        return CAMP_FILE

    return None

# *****************************************************************************
