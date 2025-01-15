# *****************************************************************************
#
# Configuration file for custom reporters.
#
# *****************************************************************************

import json

from emod_report_func import report_strain
from emod_constants import REPORTS_FILE

# *****************************************************************************


def dllcBuilder():

    # Get variables for this simulation
    # N/A

    # Dictionary to be written
    json_set = dict()

    # Custom reports object
    json_set['Custom_Reports'] = dict()

    # Configurations
    report_strain(json_set, every_timestep=True)

    #  Write file
    with open(REPORTS_FILE, 'w') as fid01:
        json.dump(json_set, fid01, sort_keys=True, indent=4)

    return None

# *****************************************************************************
