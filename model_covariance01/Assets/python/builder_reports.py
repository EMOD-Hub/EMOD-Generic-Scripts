# *****************************************************************************
#
# Configuration file for custom reporters.
#
# *****************************************************************************

from emod_report_func import build_file, write_file

# *****************************************************************************


def reportsBuilder():

    # Dictionary to be written
    json_set = build_file()

    # Configurations
    pass

    #  Write file
    write_file(json_set)

    return None

# *****************************************************************************
