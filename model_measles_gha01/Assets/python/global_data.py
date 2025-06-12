# *****************************************************************************
#
# Static data module for embedded python.
#
# *****************************************************************************

# Control params
sim_index = 0
var_params = dict()

# Filename params
demog_files = list()

demog_object = None

# Other stuff
demog_min_pop = 50

node_idval = None
demog_dict = dict()
demog_node = dict()
demog_node_map = dict()  # LGA:  [NodeIDs]
demog_rep_index = dict()  # LGA:  Output row number

first_call_bool = True
prev_proc_time = -1.0
max_node_id = 0

data_vec_time = None
data_vec_node = None
data_vec_mcw = None

adm01_list = None
nobs_vec = None

inf_dur_mean = 18.0

brate_mult_x = None
brate_mult_y = None

start_year = 2008.0
start_year_log = 2011.0
base_year = 1900.0
t_step_days = 1.0

max_clock = 180.0

# *****************************************************************************
