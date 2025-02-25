# *****************************************************************************
#
# *****************************************************************************

import os
import sys

from idmtools.core.platform_factory import Platform

# Ought to go in emodpy
LOCAL_PATH = os.path.abspath(os.path.join('..', '..', 'local_python'))
sys.path.insert(0, LOCAL_PATH)
from emod_exp import exp_from_def_file
from emod_reduce import FILENAME_ID

# *****************************************************************************

# Paths
PATH_EXP_DEF = os.path.abspath('param_dict.json')
PATH_PYTHON = os.path.abspath(os.path.join('..', 'Assets', 'python'))
PATH_DATA = os.path.abspath(os.path.join('..', 'Assets', 'data'))
PATH_EXE = os.path.abspath(os.path.join('..', '..', 'env_Debian12'))

# *****************************************************************************


# Start and experiment on COMPS
def run_sims():

    # Prepare the platform
    plat_obj = Platform(block='COMPS',
                        endpoint='https://comps.idmod.org',
                        environment='Calculon',
                        priority='Normal',
                        simulation_root='$COMPS_PATH(USER)',
                        node_group='idm_abcd',
                        num_cores='1',
                        num_retries='0',
                        exclusive='False')

    # Create experiment object
    exp_obj = exp_from_def_file(PATH_EXP_DEF, PATH_PYTHON, PATH_EXE, PATH_DATA)

    # Send experiment to COMPS; start processing
    plat_obj.run_items(exp_obj)

    # Save experiment id to file
    exp_obj.to_id_file(FILENAME_ID)
    print()
    print(exp_obj.uid.hex)

    return None

# ******************************************************************************


if (__name__ == "__main__"):

    run_sims()

# ******************************************************************************
