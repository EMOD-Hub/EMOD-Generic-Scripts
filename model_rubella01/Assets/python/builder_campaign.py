# *****************************************************************************
#
# Campaign file.
#
# *****************************************************************************

import global_data as gdata

import numpy as np

import emod_api.campaign as camp_module

from emod_camp_events import ce_br_force, ce_RI, ce_SIA
from emod_constants import CAMP_FILE

# *****************************************************************************


def campaignBuilder():

    # Variables for this simulation
    ADD_SIA = gdata.var_params['add_campaigns']
    RI_RATE = gdata.var_params['RI_rate']
    RI_XVEC = np.array(gdata.var_params['RI_rate_mult_xvals'])
    RI_YVEC = np.array(gdata.var_params['RI_rate_mult_yvals'])

    # Note: campaign module itself is the file object; no Campaign class
    ALL_NODES = gdata.demog_object.node_ids

    # Time varying birth rate
    BR_MULT_X = gdata.brate_mult_x
    BR_MULT_Y = gdata.brate_mult_y
    start_day = 365.0*(gdata.start_year-gdata.base_year)
    camp_event = ce_br_force(ALL_NODES, BR_MULT_X, BR_MULT_Y, start_day)
    camp_module.add(camp_event)

    # Routine immunization
    start_day = 365.0*(gdata.start_year-gdata.base_year+gdata.ri_offset)
    ri_yvals = np.minimum(RI_YVEC*RI_RATE, 1.0).tolist()
    ri_xvals = ((RI_XVEC-gdata.start_year-gdata.ri_offset)*365.0).tolist()
    camp_event = ce_RI(ALL_NODES, ri_xvals, ri_yvals, start_day)
    camp_module.add(camp_event)

    # SIAs
    if (ADD_SIA):
        # Catch-up
        start_day = 365.0*(gdata.start_year-gdata.base_year+gdata.ri_offset)
        camp_event = ce_SIA(ALL_NODES, start_day=start_day,
                            yrs_min=0.75, yrs_max=15.0, coverage=0.8)
        camp_module.add(camp_event)

    # End file construction
    camp_module.save(filename=CAMP_FILE)

    return None

# *****************************************************************************
