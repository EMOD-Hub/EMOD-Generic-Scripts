# *****************************************************************************
#
# Demographics file and overlays.
#
# *****************************************************************************

import os

import global_data as gdata

import numpy as np

from emod_api.demographics.Demographics import Demographics, Node

from emod_demog_func import demog_vd_calc, demog_vd_over, demog_is_over
from emod_constants import DEMOG_FILE

# *****************************************************************************


def demographicsBuilder():

    # Variables for this simulation
    PROC_DISPER = gdata.var_params['proc_overdispersion']
    SUB_ADM02 = True
    IND_RISK_VAR = gdata.var_params['ind_variance_risk']
    R0 = gdata.var_params['R0']
    LOG10_IMP = gdata.var_params['log10_import_rate']

    # Load demographic data
    dat_file = 'pop_dat_GHA.csv'
    fname_pop = os.path.join('Assets', 'data', dat_file)
    pop_input = np.loadtxt(fname_pop, dtype=int, delimiter=',')
    year_vec = pop_input[0, :] - gdata.base_year
    year_init = gdata.start_year - gdata.base_year
    pop_mat = pop_input[1:, :] + 0.1
    pop_init = [np.interp(year_init, year_vec, pop_mat[idx, :])
                for idx in range(pop_mat.shape[0])]

    # Load population data
    fname = 'demog_data_ADM02.csv'
    if (SUB_ADM02):
        fname = 'demog_data_ADM02_sub_100km.csv'
    fname = os.path.join('Assets', 'data', fname)
    with open(fname) as fid01:
        flines = [lval.strip().split(',') for lval in fid01.readlines()[1:]]
    n_nam = np.array([lval[0] for lval in flines], dtype=str)
    n_lng = np.array([lval[1] for lval in flines], dtype=float)
    n_lat = np.array([lval[2] for lval in flines], dtype=float)
    n_pop = np.array([lval[3] for lval in flines], dtype=int)
    n_pop = n_pop*(np.sum(pop_init)/np.sum(n_pop))

    # Construct node list
    list_nam = n_nam[n_pop > gdata.demog_min_pop]
    list_lat = n_lat[n_pop > gdata.demog_min_pop]
    list_lng = n_lng[n_pop > gdata.demog_min_pop]
    list_pop = n_pop[n_pop > gdata.demog_min_pop]
    node_list = [Node(lat=list_lat[k1],
                      lon=list_lng[k1],
                      pop=list_pop[k1],
                      name=list_nam[k1],
                      forced_id=(k1+1)) for k1 in range(list_nam.shape[0])]

    imp_case = np.power(10.0, LOG10_IMP)
    for node_obj in node_list:
        imp_rate = imp_case*node_obj.node_attributes.initial_population/1.0e5
        add_attrib = {'InfectivityReservoirSize': imp_rate}
        node_obj.node_attributes.extra_attributes = add_attrib

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

    # Create primary file
    ref_name = 'measles-custom'
    demog_obj = Demographics(nodes=node_list, idref=ref_name)

    # Save filename to global data for use in other functions
    gdata.demog_files.append(DEMOG_FILE)

    # Update defaults in primary file
    demog_obj.raw['Defaults']['IndividualAttributes'].clear()
    iadict = {'AcquisitionHeterogeneityVariance': IND_RISK_VAR}
    demog_obj.raw['Defaults']['IndividualAttributes'].update(iadict)

    demog_obj.raw['Defaults']['NodeAttributes'].clear()
    nadict = {'InfectivityOverdispersion': PROC_DISPER,
              'InfectivityMultiplier': 1.0}
    demog_obj.raw['Defaults']['NodeAttributes'].update(nadict)

    # Calculate vital dynamics
    vd_tup = demog_vd_calc(year_vec, year_init, pop_mat, pop_init)

    mort_year = vd_tup[0]
    mort_mat = vd_tup[1]
    age_x = vd_tup[2]
    age_y = None
    birth_rate = vd_tup[3]
    br_mult_x = vd_tup[4]
    br_mult_y = vd_tup[5]

    gdata.brate_mult_x = br_mult_x.tolist()
    gdata.brate_mult_y = br_mult_y.tolist()

    # Write vital dynamics overlay
    nfname = demog_vd_over(ref_name, node_list, birth_rate,
                           mort_year, mort_mat, age_x, age_y)
    gdata.demog_files.append(nfname)

    # Write initial susceptibility overlay
    nfname = demog_is_over(ref_name, node_list, R0, age_x, age_y)
    gdata.demog_files.append(nfname)

    # Write primary demographics file
    demog_obj.generate_file(name=DEMOG_FILE)

    # Save the demographics object for use in other functions
    gdata.demog_object = demog_obj

    return None

# *****************************************************************************
