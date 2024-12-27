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

from emod_demog_func import demog_vd_calc, demog_vd_over, demog_is_over_precalc
from emod_constants import DEMOG_FILE, MORT_XVAL

from refdat_alias import data_dict as dict_alias
from refdat_birthrate import data_dict as dict_birth
from refdat_deathrate import data_dict as dict_death

# *****************************************************************************


def demographicsBuilder():

    # Variables for this simulation
    GEOG_LIST = ['AFRO:NIGERIA']

    SUB_LGA = gdata.var_params['use_10k_res']
    START_YEAR = gdata.var_params['start_year']
    PROC_DISPER = gdata.var_params['proc_overdispersion']
    IND_RISK_VAR = gdata.var_params['ind_variance_risk']

    PATH_OVERLAY = 'demog_overlay'
    if (not os.path.exists(PATH_OVERLAY)):
        os.mkdir(PATH_OVERLAY)

    # Populate nodes in primary file
    fname = 'demog_NGA.json'
    with open(os.path.join('Assets', 'data', fname)) as fid01:
        demog_dat = json.load(fid01)

    # Aggregate population, area, lat, long, year
    list_adm02 = sorted(list(set([val.rsplit(':',1)[0] for val in demog_dat])))
    if (not SUB_LGA):
        ndemog_dat = {val: 5*[0] for val in list_adm02}
        for val in demog_dat:
            nval = val.rsplit(':',1)[0]
            ndemog_dat[nval][0] += demog_dat[val][0]
            ndemog_dat[nval][1] += demog_dat[val][1]
            ndemog_dat[nval][2] += demog_dat[val][2]*demog_dat[val][1]
            ndemog_dat[nval][3] += demog_dat[val][3]*demog_dat[val][1]
            ndemog_dat[nval][4] = demog_dat[val][4]
        demog_dat = ndemog_dat
        for val in demog_dat:
            demog_dat[val][2] /= demog_dat[val][1]
            demog_dat[val][3] /= demog_dat[val][1]

    # Add nodes
    node_id = 0
    node_list = list()
    ipop_time = dict()

    for loc_name in demog_dat:
        if (any([loc_name.startswith(val+':') or
                 val == loc_name for val in GEOG_LIST])):
            ipop_time[loc_name] = demog_dat[loc_name][4] + 0.5
            node_id = node_id + 1

            node_obj = Node(lat=demog_dat[loc_name][3],
                            lon=demog_dat[loc_name][2],
                            pop=demog_dat[loc_name][0],
                            name=loc_name,
                            forced_id=node_id,
                            area=demog_dat[loc_name][1])
            node_obj.node_attributes.infectivity_multiplier = 1.0

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
    rev_node_list = [node_obj for node_obj in node_list
                     if node_obj.node_attributes.initial_population >= gdata.demog_min_pop]
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
    demog_obj.raw['Defaults']['NodeAttributes'].clear()
    nadict = dict()
    nadict['InfectivityOverdispersion'] = PROC_DISPER
    demog_obj.raw['Defaults']['NodeAttributes'].update(nadict)

    demog_obj.raw['Defaults']['IndividualAttributes'].clear()
    iadict = dict()
    iadict['AcquisitionHeterogeneityVariance'] = IND_RISK_VAR
    demog_obj.raw['Defaults']['IndividualAttributes'].update(iadict)

    # Populate vital dynamics overlays

    # Remove mortality aliases from ref data
    ddeath_keys = list(dict_death.keys())
    for rname in ddeath_keys:
        if (rname in dict_alias):
            for rname_val in dict_alias[rname]:
                dict_death[rname_val] = dict_death[rname]
            dict_death.pop(rname)

    # Remove birth aliases from ref data
    for rname in dict_birth:
        if (rname in dict_alias):
            for rname_val in dict_alias[rname]:
                dict_birth[rname_val] = dict_birth[rname]
            dict_birth.pop(rname)

    # Create list of vital dynamics overlays
    vd_over_list = list()
    for node_dict in demog_obj.nodes:
        node_name = node_dict.name
        node_birth = None
        node_death = None

        for rname in dict_death:
            if (node_name.startswith(rname+':') or node_name == rname):
                if (node_death is None):
                    node_death = np.array(dict_death[rname])
                else:
                    raise Exception("Duplicate mortality data")

        for rname in dict_birth:
            if (node_name.startswith(rname+':') or node_name == rname):
                if (node_birth is None):
                    node_birth = np.array(dict_birth[rname])
                else:
                    raise Exception("Duplicate birthrate data")

        match_found = False
        for data_tup in vd_over_list:
            if (np.all(np.equal(data_tup[0], node_birth)) and
               np.all(np.equal(data_tup[1], node_death))):
                data_tup[2].append(node_dict)
                match_found = True
                break

        if (not match_found):
            vd_over_list.append((node_birth, node_death, [node_dict]))

    # Write vital dynamics overlays
    for k1 in range(len(vd_over_list)):
        data_tup = vd_over_list[k1]

        # Age initialization magic
        brth_rate = data_tup[0].item()
        d_rate_y = data_tup[1].tolist()
        d_rate_x = dict_death['BIN_EDGES']
        force_v = 12*[1.0]  # No seasonal forcing
        (grate, age_x, age_y) = DT._computeAgeDist(brth_rate, d_rate_x,
                                                   d_rate_y, force_v,
                                                   max_yr=100)

        # Update initial node populations; calculate target populations
        ref_nodes = [node_obj_dict.forced_id for node_obj_dict in data_tup[2]]
        for node_dict in demog_obj.nodes:
            node_name = node_dict.name
            node_id = node_dict.forced_id
            if (node_id in ref_nodes):
                ref_year = ipop_time[loc_name]
                mult_fac = grate**(START_YEAR-ref_year)
                new_pop = int(mult_fac * node_dict.node_attributes.initial_population)
                node_dict.node_attributes.initial_population = new_pop

        # Write vital dynamics overlay
        start_year_vec = [100.0]
        mort_mat = data_tup[1][..., np.newaxis]
        nfn = demog_vd_over(ref_name, data_tup[2], brth_rate, start_year_vec,
                            mort_mat, age_x, age_y, k1, d_rate_x)
        gdata.demog_files.append(nfn)

    # Load immunity mapper data
    fname = 'sus_init_NGA_{:02d}.json'.format(gdata.init_coverage)
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
