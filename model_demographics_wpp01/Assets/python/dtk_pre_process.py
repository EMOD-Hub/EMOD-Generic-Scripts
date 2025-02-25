# *****************************************************************************
#
# *****************************************************************************

from emod_preproc_func import standard_pre_process

# *****************************************************************************


def application(config_filename_in):

    # Standard pre-process function calls all other file builders
    config_filename = standard_pre_process()

    # Pre-process function needs to return config filename as a string
    return config_filename

# *****************************************************************************
