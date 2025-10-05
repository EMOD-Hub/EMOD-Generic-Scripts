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
from emod_constants import DEMOG_FILE, DEMOG_IRS

# *****************************************************************************


def demographicsBuilder():

    # Variables for this simulation
    R0 = gdata.var_params['R0']
    LOG10_IMP = gdata.var_params['log10_import_mult']
    REF_YEAR = gdata.var_params['ref_year']

    # Demographic reference data file
    dat_file = 'pop_dat_COD.csv'
    fname_pop = os.path.join('Assets', 'data', dat_file)

    # Calculate vital dynamics
    vd_tup = demog_vd_calc(fname_pop, REF_YEAR, steady_state=True)

    gdata.brate_mult_x = vd_tup[5]
    gdata.brate_mult_y = vd_tup[6]

    # Populate nodes in primary file
    node_list = list()
    node_id = 1
    imp_rate = R0/6.0 * gdata.init_pop * 1.615e-7 * np.power(10.0, LOG10_IMP)
    nname = 'EXAMPLE:A{:05d}'.format(node_id)
    node_obj = Node(lat=0.0, lon=0.0, pop=gdata.init_pop,
                    name=nname, forced_id=node_id)
    irs_dict = {DEMOG_IRS: imp_rate}
    node_obj.node_attributes.extra_attributes = irs_dict
    node_list.append(node_obj)

    # Create primary file
    ref_name = 'Demographics_Datafile'
    demog_obj = Demographics(nodes=node_list, idref=ref_name)

    # Save filename to global data for use in other functions
    gdata.demog_files.append(DEMOG_FILE)

    # Update defaults in primary file
    demog_obj.raw['Defaults']['IndividualAttributes'].clear()
    demog_obj.raw['Defaults']['NodeAttributes'].clear()

    # Write vital dynamics overlay
    nfname = demog_vd_over(ref_name, node_list, vd_tup[4],
                           vd_tup[1], vd_tup[2], vd_tup[3])
    gdata.demog_files.append(nfname)

    # Write initial susceptibility overlay
    nfname = demog_is_over(ref_name, node_list, R0, vd_tup[3])
    gdata.demog_files.append(nfname)

    # Write primary demographics file
    demog_obj.generate_file(name=DEMOG_FILE)

    # Save the demographics object for use in other functions
    gdata.demog_object = demog_obj

    return None

# *****************************************************************************
