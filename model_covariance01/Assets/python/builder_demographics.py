# *****************************************************************************
#
# Demographics file and overlays.
#
# *****************************************************************************

import global_data as gdata

from emod_api.demographics.demographics import Demographics
from emod_api.demographics.node import Node

from emod_constants import DEMOG_FILE

# *****************************************************************************


def demographicsBuilder():

    # Variables for this simulation
    IND_RISK_VAR = gdata.var_params['indiv_variance_acq']

    # Populate nodes in primary file
    node_list = list()
    node_obj = Node(lat=0.0, lon=0.0, pop=100000,
                    name='CATAN:001', forced_id=1)
    node_list.append(node_obj)

    # Create primary file
    ref_name = 'Example_Covariance_Sims'
    demog_obj = Demographics(nodes=node_list, idref=ref_name)

    # Update defaults in primary file
    demog_obj.default_node.individual_attributes.parameter_dict = dict()
    demog_obj.default_node.node_attributes.parameter_dict = dict()

    iadict = dict()
    iadict['AcquisitionHeterogeneityVariance'] = IND_RISK_VAR
    demog_obj.default_node.individual_attributes.parameter_dict.update(iadict)

    # Write primary demographics file
    demog_obj.to_file(path=DEMOG_FILE)

    # Save filename to global data for use in other functions
    gdata.demog_files.append(DEMOG_FILE)

    # Save the demographics object for use in other functions
    gdata.demog_object = demog_obj

    return None

# *****************************************************************************
