# *****************************************************************************
#
# Demographics file and overlays.
#
# *****************************************************************************

import os

import global_data as gdata

import numpy as np

from emod_api.demographics.demographics import Demographics
from emod_api.demographics.node import Node

from emod_demog_func import demog_vd_calc, demog_vd_over, demog_is_over
from emod_constants import DEMOG_FILE, DEMOG_IRS

# *****************************************************************************


def demographicsBuilder():

    # Variables for this simulation
    R0 = gdata.var_params['R0']
    LOG10_IMP = gdata.var_params['log10_import_mult']
    SS_DEMOG = gdata.var_params['steady_state_demog']
    REF_CODE = gdata.var_params['demog_set']
    NUM_NODES = gdata.var_params['num_nodes']

    # Demographic reference data file
    dat_file = f'pop_dat_{REF_CODE}.csv'
    fname_pop = os.path.join('Assets', 'data', dat_file)

    # Calculate vital dynamics
    vd_tup = demog_vd_calc(fname_pop, gdata.start_year, steady_state=SS_DEMOG)

    gdata.brate_mult_x = vd_tup[5]
    gdata.brate_mult_y = vd_tup[6]

    # Populate nodes in primary file
    node_list = list()
    pop_nodes = int(gdata.init_pop/NUM_NODES)
    imp_rate = pop_nodes * 1.615e-7 * np.power(10.0, LOG10_IMP)
    for nid in range(NUM_NODES):
        nname = REF_CODE + '_{:04d}'.format(nid+1)
        node_obj = Node(lat=0.0, lon=0.0, pop=pop_nodes,
                        name=nname, forced_id=(nid+1))
        irs_dict = {DEMOG_IRS: imp_rate}
        node_obj.node_attributes.extra_attributes = irs_dict
        node_list.append(node_obj)

    # Create primary file
    ref_name = 'Demographics_Datafile'
    demog_obj = Demographics(nodes=node_list, idref=ref_name,
                             set_defaults=False)

    # Update defaults in primary file
    demog_obj.default_node.individual_attributes.parameter_dict = dict()
    demog_obj.default_node.node_attributes.parameter_dict = dict()

    # Write primary demographics file
    demog_obj.to_file(path=DEMOG_FILE)

    # Save filename to global data for use in other functions
    gdata.demog_files.append(DEMOG_FILE)

    # Write initial susceptibility overlay
    nfname = demog_is_over(ref_name, node_list, R0, vd_tup[3])
    gdata.demog_files.append(nfname)

    # Write vital dynamics overlay
    nfname = demog_vd_over(ref_name, node_list, vd_tup[4],
                           vd_tup[1], vd_tup[2], vd_tup[3])
    gdata.demog_files.append(nfname)

    # Save the demographics object for use in other functions
    gdata.demog_object = demog_obj

    return None

# *****************************************************************************
