# *****************************************************************************
#
# Configuration file for custom reporters.
#
# *****************************************************************************

from emod_report_func import build_file, write_file, report_strain

import global_data as gdata

# *****************************************************************************


def reportsBuilder():

    # Dictionary to be written
    json_set = build_file()

    # Configurations
    report_strain(json_set, time_start=gdata.start_log)

    #  Write file
    write_file(json_set)

    return None

# *****************************************************************************
