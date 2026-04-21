# *****************************************************************************
#
# Configuration file for custom reporters.
#
# *****************************************************************************

from emod_report_func import build_file, write_file, report_strain

# *****************************************************************************


def reportsBuilder():

    # Dictionary to be written
    json_set = build_file()

    # Configurations
    report_strain(json_set, every_timestep=True)

    #  Write file
    write_file(json_set)

    return None

# *****************************************************************************
