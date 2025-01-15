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
inproc_abort_time = 43615
inproc_abort_bool = True

demog_node = dict()
demog_min_pop = 50
demog_node_map = dict() # LGA: [NodeIDs]
demog_rep_index = dict() # LGA: Output row number

brate_mult_x = None
brate_mult_y = None

init_ob_thresh = 10000

seed_inf_dt = 60.0
seed_inf_num = 25.0

inf_dur_mean = 24.0
inf_dur_std = 11.3

boxes_nopv2 = 2
boxes_sabin2 = 7
rev_nopv2 = 0.0
rev_sabin2 = 0.0
init_coverage = 50

base_year = 1900.0
t_step_days = 5.0

max_clock = 180.0

# *****************************************************************************
