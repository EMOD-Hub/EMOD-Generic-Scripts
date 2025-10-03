# *****************************************************************************
#
# Demographics file and overlays.
#
# *****************************************************************************

import os

import global_data as gdata

import numpy as np

from emod_api.demographics.Demographics import Demographics, Node

from emod_demog_func import demog_vd_calc, demog_vd_over
from emod_constants import DEMOG_FILE, BASE_YEAR

# *****************************************************************************


def demographicsBuilder():

    # Variables for this simulation
    START_YEAR = gdata.var_params['start_year']
    POP_DAT_STR = gdata.var_params['pop_dat_file']

    # Demographic reference data file
    dat_file = 'pop_dat_{:s}.csv'.format(POP_DAT_STR)
    fname_pop = os.path.join('Assets', 'data', dat_file)

    # Calculate vital dynamics
    vd_tup = demog_vd_calc(fname_pop, START_YEAR)

    gdata.init_pop = vd_tup[0]
    gdata.brate_mult_x = vd_tup[5]
    gdata.brate_mult_y = vd_tup[6]

    # Populate nodes in primary file
    node_list = list()
    node_obj = Node(lat=0.0, lon=0.0, pop=gdata.init_pop,
                    name=POP_DAT_STR, forced_id=1)
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

    # Write primary demographics file
    demog_obj.generate_file(name=DEMOG_FILE)

    # Save the demographics object for use in other functions
    gdata.demog_object = demog_obj

    return None

# *****************************************************************************
