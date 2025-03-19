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

from emod_demog_func import demog_vd_calc, demog_vd_over, demog_ah_over, \
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

    # Populate nodes in primary file
    fname = 'demog_data.json'
    with open(os.path.join('Assets', 'data', fname)) as fid01:
        demog_dat = json.load(fid01)

    # Reference population structure
    pop_ratio = dict()
    vd_data = dict()
    for adm00 in gdata.targ_adm00:
        cname = adm00.split(':')[-1]
        dat_file = 'pop_dat_{:s}.csv'.format(cname)
        fname_pop = os.path.join('Assets', 'data', dat_file)
        pop_input = np.loadtxt(fname_pop, dtype=int, delimiter=',')

        year_vec = pop_input[0, :] - gdata.base_year
        year_init = START_YEAR - gdata.base_year
        pop_mat = pop_input[1:, :] + 0.1
        pop_init = [np.interp(year_init, year_vec, pop_mat[idx, :])
                    for idx in range(pop_mat.shape[0])]

        # Vital dynamics data
        vd_data[adm00] = demog_vd_calc(year_vec, year_init, pop_mat, pop_init)

        # Population fudge factor
        file_pop = sum([demog_dat[val][0] for val in demog_dat
                        if val.startswith(adm00+':') or val == adm00])
        ref_pop = np.sum(pop_init)
        pop_ratio[cname] = ref_pop/file_pop


    # Aggregate population, area, lat, long
    stem01 = [val.rsplit(':', 1)[0] for val in demog_dat]
    list_adm02 = sorted(list(set(stem01)))
    if (not SUB_ADM02):
        ndemog_dat = {val: 4*[0] for val in list_adm02}
        for val in demog_dat:
            nval = val.rsplit(':', 1)[0]
            ndemog_dat[nval][0] += demog_dat[val][0]
            ndemog_dat[nval][1] += demog_dat[val][1]
            ndemog_dat[nval][2] += demog_dat[val][2]*demog_dat[val][1]
            ndemog_dat[nval][3] += demog_dat[val][3]*demog_dat[val][1]
        demog_dat = ndemog_dat
        for val in demog_dat:
            demog_dat[val][2] /= demog_dat[val][1]
            demog_dat[val][3] /= demog_dat[val][1]

    # Add nodes
    node_id = 0
    node_list = list()

    for loc_name in demog_dat:
        if (any([loc_name.startswith(val+':') or
                 val == loc_name for val in gdata.targ_adm00])):
            node_id = node_id + 1
            cname = loc_name.split(':')[1]
            pop_ff = pop_ratio[cname]
            node_obj = Node(lat=demog_dat[loc_name][3],
                            lon=demog_dat[loc_name][2],
                            pop=demog_dat[loc_name][0]*pop_ff,
                            name=loc_name,
                            forced_id=node_id,
                            area=demog_dat[loc_name][1])
            node_list.append(node_obj)

    # Node name bookkeeping
    nname_dict = {node_obj.name: node_obj.forced_id for node_obj in node_list}
    rep_groups = {nrep: [nname_dict[val] for val in nname_dict.keys()
                         if val.startswith(nrep+':') or val == nrep]
                  for nrep in list_adm02}
    gdata.demog_node_map = rep_groups

    node_rep_dict = {val[0]: val[1] for val in
                     zip(list_adm02, range(len(list_adm02)))}
    gdata.demog_rep_index = node_rep_dict

    # Prune small nodes
    rev_node_list = list()
    for node_obj in node_list:
        n_init_pop = node_obj.node_attributes.initial_population
        if (n_init_pop >= gdata.demog_min_pop):
            rev_node_list.append(node_obj)
    node_list = rev_node_list

    node_name_dict = {node_obj.name: node_obj.forced_id
                      for node_obj in node_list}
    gdata.demog_node = node_name_dict

    # Create primary file
    ref_name = 'polio-custom'
    demog_obj = Demographics(nodes=node_list, idref=ref_name)

    # Save filename to global data for use in other functions
    gdata.demog_files.append(DEMOG_FILE)

    # Update defaults in primary file
    demog_obj.raw['Defaults']['IndividualAttributes'].clear()
    demog_obj.raw['Defaults']['NodeAttributes'].clear()

    nadict = dict()
    nadict['InfectivityOverdispersion'] = PROC_DISPER
    demog_obj.raw['Defaults']['NodeAttributes'].update(nadict)

    # Spatial R0 variability
    fname = os.path.join('Assets', 'data', 'cgf_index_underweight.json')
    with open(fname) as fid01:
        dict_ahv = json.load(fid01)

    k1 = 0
    for reg_name in dict_ahv:
        p_val = dict_ahv[reg_name]
        ah_val = AH_VAR_SCALE
        ahv_mean = 0.21172
        r0_val = 1.0/(1.0+np.exp(R0_SCALE*(ahv_mean-p_val))) + R0_MIN_M

        n_list = [n_obj for n_obj in node_list
                  if (n_obj.name == reg_name or
                      n_obj.name.startswith(reg_name+':'))]
        nfname = demog_ah_over(ref_name, n_list, r0_val, ah_val, k1)
        gdata.demog_files.append(nfname)
        k1 = k1 + 1

    # Write vital dynamics overlay
    k1 = 0
    for adm00 in vd_data:
        vd_tup = vd_data[adm00]
        mort_year = vd_tup[0]
        mort_mat = vd_tup[1]
        age_x = vd_tup[2]
        age_y = None
        birth_rate = vd_tup[3]

        n_list = [n_obj for n_obj in node_list
                  if (n_obj.name == adm00 or
                      n_obj.name.startswith(adm00+':'))]

        nfname = demog_vd_over(ref_name, n_list, birth_rate,
                               mort_year, mort_mat, age_x, age_y, k1)
        gdata.demog_files.append(nfname)
        k1 = k1 + 1

        node_ids = [node_obj.forced_id for node_obj in n_list]
        br_tup = ((vd_tup[4]).tolist(), (vd_tup[5]).tolist(), node_ids)
        gdata.brate_mult_tup.append(br_tup)

    ## Birth rate indexing
    #fname = 'cbr_NGA.json'
    #with open(os.path.join('Assets', 'data', fname)) as fid01:
    #    cbr_mult_dict = json.load(fid01)

    # Write vital dynamics overlay

    #for reg_name in cbr_mult_dict:
    #    n_cbr = cbr_mult_dict[reg_name]*birth_rate
    #    mort_vec = np.array([np.interp(year_init, mort_year, mort_mat[idx, :])
    #                         for idx in range(mort_mat.shape[0])])
    #    mort_vec = mort_vec.tolist()
    #    f_vec = 12*[1.0]  # No seasonal forcing
    #    (_, age_x, age_y) = DT._computeAgeDist(n_cbr, MORT_XVAL,
    #                                           mort_vec, f_vec)

    #    n_list = [n_obj for n_obj in node_list
    #              if (n_obj.name == reg_name or
    #                  n_obj.name.startswith(reg_name+':'))]

    # Load immunity mapper data
    fname = 'sus_init_{:02d}.json'.format(gdata.init_coverage)
    with open(os.path.join('Assets', 'data', fname)) as fid01:
        isus_dat = json.load(fid01)

    isus_time = np.array(isus_dat['time'])
    isus_name = np.array(isus_dat['name'])
    isus_ages = np.array(isus_dat['ages'])
    isus_data = np.array(isus_dat['data'])

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
