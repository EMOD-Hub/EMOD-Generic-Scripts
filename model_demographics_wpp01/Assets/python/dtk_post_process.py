# *****************************************************************************
#
# *****************************************************************************

import json

import global_data as gdata

from emod_postproc_func import post_proc_poppyr
from emod_constants import O_FILE

# *****************************************************************************


def application(output_path):

    # Prep output dictionary
    SIM_IDX = gdata.sim_index
    key_str = '{:05d}'.format(SIM_IDX)
    parsed_dat = {key_str: dict()}

    # Sample population pyramid every year
    post_proc_poppyr(output_path, parsed_dat[key_str])

    # Write output dictionary
    with open(O_FILE, 'w') as fid01:
        json.dump(parsed_dat, fid01)

    return None

# *****************************************************************************
