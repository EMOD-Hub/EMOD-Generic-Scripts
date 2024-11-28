# *****************************************************************************
#
# *****************************************************************************

import json
import os
import sys

from idmtools.core.id_file import read_id_file

# Ought to go in emodpy
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'local_python')))
from emod_reduce import get_sim_files, pool_manager
from py_assets_common.emod_constants import COMPS_ID_FILE, D_FILE

# *****************************************************************************


def get_data_brick():

    # Get Experiment ID
    (exp_id, _, _, _) = read_id_file(COMPS_ID_FILE)

    # Reduce output and write data brick
    pool_manager(exp_id)
    #data_brick = get_sim_files(exp_id)
    #with open(D_FILE, 'w') as fid01:
        #json.dump(data_brick, fid01)

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    get_data_brick()

# *****************************************************************************
