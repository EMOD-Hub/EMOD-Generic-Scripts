# *****************************************************************************
#
# *****************************************************************************

import json
import os
import sys

from idmtools.core.id_file import read_id_file

# Ought to go in emodpy
LOCAL_PATH = os.path.abspath(os.path.join('..', '..', 'local_python'))
sys.path.insert(0, LOCAL_PATH)
from emod_reduce import get_sim_files

# *****************************************************************************


def get_data_brick():

    # Get Experiment ID
    (exp_id, _, _, _) = read_id_file('COMPS_ID.id')

    # Reduce output and write data brick
    data_brick = get_sim_files(exp_id, LOCAL_PATH)
    with open('data_brick.json', 'w') as fid01:
        json.dump(data_brick, fid01)

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    get_data_brick()

# *****************************************************************************
