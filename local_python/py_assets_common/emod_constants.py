# *****************************************************************************
#
# *****************************************************************************

API_MIN = '2.0.0'

COMPS_ID_FILE = 'COMPS_ID.id'
COMPS_URL = 'https://comps.idmod.org'
COMPS_SU_EXE = 'Execution'
COMPS_SU_ENV = 'Environment'
COMPS_SU_CMD = 'export PYTHONPATH="$PYTHONPATH:$COMPS_WORKING_DIR/Assets";' + \
               '/py_env/bin/python3 /mnt/idm2/comps/scripts/startup.py'
COMPS_SU_DICT = {'COMPS_STARTUP_CMD': COMPS_SU_CMD}

LOCAL_EXP_ROOT = 'docker_suites'
LOCAL_EXP_DIR = 'LOCAL_EXP_PATH.txt'

PY_PATH = r'python3'

DOCK_PACK = r'mesokurtic/emodpy:2.1.0.03'

VE_PY_PATHS = ['/py_env/lib/python3.9/site-packages/',
               '/py_env/lib/python3.10/site-packages/',
               '/py_env/lib/python3.11/site-packages/',
               '/py_env/lib/python3.12/site-packages/',
               '/py_env/lib/python3.13/site-packages/']

ID_OS = 'EMOD_OS.id'
ID_EXE = 'EMOD_EXE.id'
ID_ENV = 'EMOD_ENV.id'
ID_SCHEMA = 'EMOD_SCHEMA.id'

P_FILE = 'param_dict.json'
D_FILE = 'data_brick.json'
I_FILE = 'idx_str_file.txt'
C_FILE = 'config.json'
O_FILE = 'parsed_out.json'

EXP_C = 'EXP_CONSTANT'
EXP_V = 'EXP_VARIABLE'
EXP_O = 'EXP_OPTIMIZE'
EXP_NAME = 'EXP_NAME'
NUM_SIMS = 'NUM_SIMS'
NUM_ITER = 'NUM_ITER'

SPATH = 'Assets/schema.json'

CAMP_FILE = 'campaign.json'
REPORTS_FILE = 'reports.json'
DEMOG_FILE = 'demographics.json'
RST_FILE = 'ReportStrainTracking.csv'
SQL_FILE = 'simulation_events.db'

REPORTS_KEY = 'Reports'

PATH_OVERLAY = 'demog_overlay'

DIST_CMPLX = 'DISTRIBUTION_COMPLEX'
DIST_GAMMA = 'GAMMA_DISTRIBUTION'
DIST_GAUSS = 'GAUSSIAN_DISTRIBUTION'

SQL_TIME = 0
SQL_NODE = 2
SQL_AGENT = 3
SQL_MCW = 4
SQL_AGE = 5
SQL_LABEL = 6

RST_TIME = 0
RST_NODE = 1
RST_CLADE = 2
RST_GENOME = 3
RST_TOT_INF = 4
RST_CONT_INF = 5
RST_CONT_TOT = 6
RST_NEW_INF = 7

YR_DAYS = 365
MO_DAYS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

CBR_VEC = 'cbr_vec'
R0_VEC = 'r0_vec'
R0_TIME = 'r0_time'
POP_PYR = 'pyr_data'
INF_FRAC = 'inf_frac'
CAMP_COST = 'camp_cost'
NODE_POP_STR = 'node_pop'
NODE_IDS_STR = 'node_ids'

MORT_XVAL = [0.6, 1829.5, 1829.6, 3659.5, 3659.6, 5489.5,
             5489.6, 7289.5, 7289.6, 9119.5, 9119.6, 10949.5,
             10949.6, 12779.5, 12779.6, 14609.5, 14609.6, 16439.5,
             16439.6, 18239.5, 18239.6, 20069.5, 20069.6, 21899.5,
             21899.6, 23729.5, 23729.6, 25559.5, 25559.6, 27389.5,
             27389.6, 29189.5, 29189.6, 31019.5, 31019.6, 32849.5,
             32849.6, 34679.5, 34679.6, 36509.5, 36509.6, 38339.5]

MAX_DAILY_MORT = 0.01

POP_AGE_DAYS = [0, 1825, 3650, 5475, 7300, 9125, 10950, 12775,
                14600, 16425, 18250, 20075, 21900, 23725, 25550,
                27375, 29200, 31025, 32850, 34675, 36500]

AGE_KEY_LIST = ['<5', '5-9', '10-14', '15-19', '20-24', '25-29', '30-34',
                '35-39', '40-44', '45-49', '50-54', '55-59', '60-64', '65-69',
                '70-74', '75-79', '80-84', '85-89', '90-94', '95-99']

CLR_M = [70/255, 130/255, 180/255]
CLR_F = [238/255, 121/255, 137/255]

# *****************************************************************************
