# *****************************************************************************
#
# Campaign file.
#
# *****************************************************************************

import global_data as gdata

import emod_api.campaign as camp_module

from emod_camp_events import ce_br_force, ce_RI, ce_SIA
from emod_constants import CAMP_FILE, BASE_YEAR

# *****************************************************************************


def campaignBuilder():

    # Variables for this simulation
    ADD_SIA = gdata.var_params['add_campaigns']
    RI_RATE = gdata.var_params['RI_rate']

    # Note: campaign module itself is the file object; no Campaign class
    ALL_NODES = gdata.demog_object.node_ids

    # Time varying birth rate
    BR_MULT_X = gdata.brate_mult_x
    BR_MULT_Y = gdata.brate_mult_y
    start_day = 365.0*(gdata.start_year-BASE_YEAR)
    camp_event = ce_br_force(ALL_NODES, BR_MULT_X, BR_MULT_Y, start_day)
    camp_module.add(camp_event)

    # Routine immunization
    start_day = 365.0*(gdata.start_year-BASE_YEAR+gdata.ri_offset)
    camp_event = ce_RI(ALL_NODES, coverage=RI_RATE, start_day=start_day)
    camp_module.add(camp_event)

    # SIAs
    if (ADD_SIA):
        # Catch-up
        start_day = 365.0*(gdata.start_year-BASE_YEAR+gdata.ri_offset)
        camp_event = ce_SIA(ALL_NODES, coverage=0.8, start_day=start_day,
                            yrs_min=0.75, yrs_max=15.0)
        camp_module.add(camp_event)

    # End file construction
    camp_module.save(filename=CAMP_FILE)

    return None

# *****************************************************************************
