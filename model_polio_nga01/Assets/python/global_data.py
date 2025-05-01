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
inproc_dt = 90

demog_min_pop = 75

node_idval = dict()
adm02_idlist  = dict()
adm01_idlist  = dict()

brate_mult_tup_list = list()

init_ob_thresh = 10000

targ_adm00 = ['AFRO:CAMEROON',
              'AFRO:CHAD',
              'AFRO:NIGER',
              'AFRO:NIGERIA']

seed_sets = [('AFRO:NIGERIA:JIGAWA:BIRINIWA', 2017.9),
             ('AFRO:NIGERIA:ZAMFARA:SHINKAFI', 2020.9)]

r0_mid_val = 0.20  #0.21172

nopv2_sia_take_fac = 0.7

seed_inf_dt = 60.0
seed_inf_num = 25.0

inf_tot_shape = 1.0
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
