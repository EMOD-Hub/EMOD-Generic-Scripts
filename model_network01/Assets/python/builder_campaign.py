# *****************************************************************************
#
# Campaign file.
#
# *****************************************************************************

import emod_api.campaign as camp_module

from emod_camp_events import ce_import_pressure
from emod_constants import CAMP_FILE

# *****************************************************************************


def campaignBuilder():

    # Variables for this simulation
    # N/A

    # Note: campaign module itself is the file object; no Campaign class

    # Import pressure
    camp_event = ce_import_pressure([1], duration=5.0)
    camp_module.add(camp_event)

    # End file construction
    camp_module.save(filename=CAMP_FILE)

    return None

# *****************************************************************************
