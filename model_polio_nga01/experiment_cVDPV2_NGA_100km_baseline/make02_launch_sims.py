# *****************************************************************************
#
# *****************************************************************************

import os
import sys

# Ought to go in emodpy
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'local_python')))
from emod_exp import start_exp
from py_assets_common.emod_constants import P_FILE

# *****************************************************************************

# Paths
PATH_EXP_DEF = os.path.abspath(P_FILE)
PATH_PYTHON = os.path.abspath(os.path.join('..', 'Assets', 'python'))
PATH_DATA = os.path.abspath(os.path.join('..', 'Assets', 'data'))

# *****************************************************************************


if (__name__ == "__main__"):

    start_exp(PATH_PYTHON, PATH_DATA, PATH_EXP_DEF, num_cores=2)

# ******************************************************************************
