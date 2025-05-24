# *****************************************************************************
#
# *****************************************************************************

import os
import sys

from idmtools.core.id_file import read_id_file

# Ought to go in emodpy
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'local_python')))
from emod_reduce import get_data_brick, pool_manager
from py_assets_common.emod_constants import COMPS_ID_FILE

# *****************************************************************************


if (__name__ == "__main__"):

    # get_data_brick()
    (exp_id, _, _, _) = read_id_file(COMPS_ID_FILE)
    pool_manager(exp_id)

# *****************************************************************************
