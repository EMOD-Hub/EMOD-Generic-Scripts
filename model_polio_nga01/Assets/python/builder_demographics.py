# *****************************************************************************
#
# Demographics file and overlays.
#
# *****************************************************************************

import os
import json

import global_data as gdata

import numpy as np

from emod_api.demographics.Demographics import Demographics, Node
from emod_api.demographics import DemographicsTemplates as DT

from emod_demog_func import demog_vd_calc, demog_vd_over, demog_r0mult_over, \
                            demog_is_over_precalc
from emod_constants import DEMOG_FILE, MORT_XVAL

# *****************************************************************************


def demographicsBuilder():

    # Variables for this simulation
    SUB_ADM02 = gdata.var_params['use_10k_res']

    START_YEAR = gdata.var_params['start_year']
    PROC_DISPER = gdata.var_params['proc_overdispersion']
    AH_VAR_SCALE = gdata.var_params['ind_variance_risk']
    R0_SCALE = gdata.var_params['R0_sig_scale']
    R0_MIN_M = gdata.var_params['R0_min_mult']

    # Load ISO3 codes dict
    fname = os.path.join('Assets', 'data', 'dict_iso3.json')
    with open(fname) as fid01:
        iso3_dict = json.load(fid01)

    # Load population data
    fname = 'demog_data_POLIS_ADM02.csv'
    if (SUB_ADM02):
        fname = 'demog_data_POLIS_ADM02_sub_100km.csv'
    fname = os.path.join('Assets', 'data', fname)
    with open(fname) as fid01:
        flines = [lval.strip().split(',') for lval in fid01.readlines()[1:]]
    n_nam = np.array([lval[0] for lval in flines], dtype=str)
    n_lng = np.array([lval[1] for lval in flines], dtype=float)
    n_lat = np.array([lval[2] for lval in flines], dtype=float)
    n_pop = np.array([lval[3] for lval in flines], dtype=int)
    p_mul = np.zeros(n_pop.shape[0], dtype=float)

    # Load demographic data
    vd_data = dict()
    for cname in gdata.targ_adm00:
        dat_file = 'pop_dat_{:s}.csv'.format(iso3_dict[cname])
        fname = os.path.join('Assets', 'data', dat_file)
        pop_input = np.loadtxt(fname, dtype=int, delimiter=',')

        year_vec = pop_input[0, :] - gdata.base_year
        year_init = START_YEAR - gdata.base_year
        pop_mat = pop_input[1:, :] + 0.1
        pop_init = [np.interp(year_init, year_vec, pop_mat[idx, :])
                    for idx in range(pop_mat.shape[0])]

        n_idx = np.array([n_val.startswith(cname+':') for n_val in n_nam])
        pol_pop = np.sum(n_pop[n_idx])
        ref_pop = np.sum(pop_init)
        p_mul[n_idx] = ref_pop/pol_pop

        vd_data[cname] = demog_vd_calc(year_vec, year_init, pop_mat, pop_init)

    # Construct node list
    p_mod = n_pop*p_mul
    list_nam = n_nam[p_mod > gdata.demog_min_pop]
    list_lat = n_lat[p_mod > gdata.demog_min_pop]
    list_lng = n_lng[p_mod > gdata.demog_min_pop]
    list_pop = p_mod[p_mod > gdata.demog_min_pop]
    node_list = [Node(lat=list_lat[k1],
                      lon=list_lng[k1],
                      pop=list_pop[k1],
                      name=list_nam[k1],
                      forced_id=(k1+1)) for k1 in range(list_nam.shape[0])]

    # Nameset for admin02
    adm02_subset = list_nam.tolist()
    if (SUB_ADM02):
        adm02_subset = [n_val.rsplit(':', 1)[0] for n_val in list_nam]
    list_adm02 = sorted(list(set(adm02_subset)))
    adm01_subset = [n_val.rsplit(':', 1)[0] for n_val in list_adm02]
    list_adm01 = sorted(list(set(adm01_subset)))

    # Node name bookkeeping
    node_dict = {nobj.name: nobj.forced_id for nobj in node_list}
    gdata.node_idval = node_dict

    adm02_dict = {adm02: [node_dict[nname] for nname in node_dict.keys()
                          if nname.startswith(adm02+':') or nname == adm02]
                  for adm02 in list_adm02}
    gdata.adm02_idlist = adm02_dict

    adm01_dict = {adm01: [node_dict[nname] for nname in node_dict.keys()
                          if nname.startswith(adm01+':') or nname == adm01]
                  for adm01 in list_adm01}
    gdata.adm01_idlist = adm01_dict

    node_rep_dict = {val[0]: val[1] for val in
                     zip(list_adm02, range(len(list_adm02)))}
    gdata.demog_rep_index = node_rep_dict

    # Create primary file
    ref_name = 'polio-custom'
    demog_obj = Demographics(nodes=node_list, idref=ref_name)

    # Save filename to global data for use in other functions
    gdata.demog_files.append(DEMOG_FILE)

    # Update defaults in primary file
    demog_obj.raw['Defaults']['IndividualAttributes'].clear()
    iadict = {'AcquisitionHeterogeneityVariance': AH_VAR_SCALE}
    demog_obj.raw['Defaults']['IndividualAttributes'].update(iadict)

    demog_obj.raw['Defaults']['NodeAttributes'].clear()
    nadict = {'InfectivityOverdispersion': PROC_DISPER}
    demog_obj.raw['Defaults']['NodeAttributes'].update(nadict)

    # Spatial R0 variability
    fname = os.path.join('Assets', 'data', 'cgf_index_underweight.json')
    with open(fname) as fid01:
        dict_r0_effect = json.load(fid01)

    k1 = 0
    for reg_name in dict_r0_effect:
        p_val = dict_r0_effect[reg_name]
        r0_val = 1.0/(1.0+np.exp(R0_SCALE*(gdata.r0_mid_val-p_val))) + R0_MIN_M

        n_list = [n_obj for n_obj in node_list
                  if (n_obj.name.startswith(reg_name+':'))]

        if (n_list):
            nfname = demog_r0mult_over(ref_name, n_list, r0_val, k1)
            gdata.demog_files.append(nfname)
            k1 = k1 + 1

    # Birth rate indexing
    fname = 'demog_data_cbr_mult.json'
    with open(os.path.join('Assets', 'data', fname)) as fid01:
        cbr_mult_dict = json.load(fid01)

    # Write vital dynamics overlay
    k1 = 0
    for cname in vd_data:
        vd_tup = vd_data[cname]
        mort_year = vd_tup[0]
        mort_mat = vd_tup[1]
        age_x = vd_tup[2]
        age_y = None
        birth_rate = vd_tup[3]

        node_ids = [n_obj.forced_id for n_obj in node_list
                    if (n_obj.name.startswith(cname+':'))]
        br_tup = ((vd_tup[4]).tolist(), (vd_tup[5]).tolist(), node_ids)
        gdata.brate_mult_tup_list.append(br_tup)

        for reg_name in cbr_mult_dict:
            if (not reg_name.startswith(cname+':')):
                continue

            n_list = [n_obj for n_obj in node_list
                      if (n_obj.name.startswith(reg_name+':'))]
            if (n_list):
                n_cbr = cbr_mult_dict[reg_name]*birth_rate
                year_init = START_YEAR - gdata.base_year
                f_vec = 12*[1.0]  # No seasonal forcing
                mort_vec = [np.interp(year_init, mort_year, mort_mat[idx, :])
                            for idx in range(mort_mat.shape[0])]

                (_, age_x, age_y) = DT._computeAgeDist(n_cbr, MORT_XVAL,
                                                       mort_vec, f_vec)

                nfname = demog_vd_over(ref_name, n_list, n_cbr,
                                       mort_year, mort_mat, age_x, age_y, k1)
                gdata.demog_files.append(nfname)
                k1 = k1 + 1

    # Load immunity mapper data
    fname = 'sus_init_{:02d}.json'.format(gdata.init_coverage)
    with open(os.path.join('Assets', 'data', fname)) as fid01:
        isus_dat = json.load(fid01)

    isus_time = np.array(isus_dat['time'])
    isus_name = np.array(isus_dat['name'])
    isus_ages = np.array(isus_dat['ages'])
    isus_data = np.array(isus_dat['data'])

    use_reg = np.array([any([n_val.startswith(adm00+':')
                             for adm00 in gdata.targ_adm00])
                        for n_val in isus_name], dtype=bool)
    isus_name = isus_name[use_reg]
    isus_data = isus_data[use_reg, :, :]

    # Create list of initial susceptibility overlays
    is_over_list = list()
    start_time = 365.0*(START_YEAR-gdata.base_year)
    for node_dict in demog_obj.nodes:
        node_name = node_dict.name
        node_initsus = None

        for k1 in range(isus_name.shape[0]):
            rname = isus_name[k1]
            if (node_name.startswith(rname+':') or node_name == rname):
                if (node_initsus is None):
                    node_initsus = list()
                    for k2 in range(isus_ages.shape[0]):
                        nis_dat = isus_data[k1, k2, :]
                        ipdat = np.interp(start_time, isus_time, nis_dat)
                        node_initsus.append(ipdat)
                    node_initsus = np.array(node_initsus)
                else:
                    raise Exception("Duplicate susceptibility data")

        match_found = False
        for data_tup in is_over_list:
            if (np.all(np.equal(data_tup[0], node_initsus))):
                data_tup[1].append(node_dict)
                match_found = True
                break

        if (not match_found):
            is_over_list.append((node_initsus, [node_dict]))

    # Write susceptibility overlays
    for k1 in range(len(is_over_list)):
        data_tup = is_over_list[k1]
        isus_x = [0.0] + (365.0*(isus_ages+3.0)/12.0).tolist() + [365.0*10.0]
        isus_y = [1.0] + data_tup[0].tolist() + [0.0]

        nfname = demog_is_over_precalc(ref_name, data_tup[1],
                                       isus_x, isus_y, k1)
        gdata.demog_files.append(nfname)

    # Write primary demographics file
    demog_obj.generate_file(name=DEMOG_FILE)

    # Save the demographics object for use in other functions
    gdata.demog_object = demog_obj

    return None

# *****************************************************************************
