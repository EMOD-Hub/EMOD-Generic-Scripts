# *****************************************************************************
#
# Configuration file for simulation.
#
# *****************************************************************************

import numpy as np

import global_data as gdata

from emod_constants import CAMP_FILE, REPORTS_FILE

# *****************************************************************************

def max_coeff_ref(exp_vals):

    if (np.min(exp_vals)<0.0 or np.max(exp_vals)>8.0):
        raise Exception('Network exponent out of range.')

    x_ref = np.array([0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 7, 8])
    y_ref = np.array([-2.794, -1.298, 0.155, 1.528, 2.797, 6.924,
                      9.774, 12.22, 14.44, 16.56, 18.65, 20.75])

    max_coeffs = np.interp(exp_vals, x_ref, y_ref).tolist()

    return max_coeffs

# *****************************************************************************


def update_config_obj(config):

    # Variables for this simulation
    R0 = gdata.var_params['R0']
    R0_OPV = gdata.var_params['R0_OPV']
    NOPV_R0_MULT = gdata.var_params['R0_nOPV_mult']

    BI_STD = gdata.var_params['base_inf_stddev_mult']

    ID_MEAN = gdata.var_params['inf_duration_mean']
    ID_STD = gdata.var_params['inf_dur_stddev_mult']

    NI_LN_MULT = gdata.var_params['net_inf_ln_mult']
    NI_POWER = gdata.var_params['net_inf_power']
    NI_MAXFRAC = gdata.var_params['net_inf_maxfrac']

    OPV_REV = gdata.var_params['OPV_rev_prob']
    OPV_BOXES = gdata.var_params['OPV_compartments']

    NOPV_REV = gdata.var_params['nOPV_rev_prob']
    NOPV_BOXES = gdata.var_params['nOPV_compartments']

    RUN_NUM = gdata.var_params['run_number']

    TIME_START = gdata.var_params['start_time']
    TIME_DELTA = gdata.var_params['num_tsteps']

    AGENT_RATE = gdata.var_params['agent_rate']

    CORR_ACQ_TRANS = gdata.var_params['corr_acq_trans']

    MAX_CLOCK = gdata.var_params['max_clock_minutes']

    LABEL_MUTES = gdata.var_params['label_by_mutator']

    # Config parameters object (read only dictionary)
    cp = config.parameters

    # Random number seed
    cp.Run_Number = RUN_NUM

    # Time
    cp.Start_Time = gdata.start_off + TIME_START
    cp.Simulation_Duration = TIME_DELTA

    cp.Enable_Termination_On_Total_Wall_Time = 1
    cp.Wall_Time_Maximum_In_Minutes = MAX_CLOCK

    # Intrahost
    if(BI_STD > 0.01):
        cp.Base_Infectivity_Distribution = 'GAMMA_DISTRIBUTION'
        cp.Base_Infectivity_Scale = R0/ID_MEAN*BI_STD*BI_STD
        cp.Base_Infectivity_Shape = 1.0/BI_STD/BI_STD
    else:
        cp.Base_Infectivity_Distribution = 'CONSTANT_DISTRIBUTION'
        cp.Base_Infectivity_Constant = R0/ID_MEAN

    cp.Incubation_Period_Distribution = 'GAUSSIAN_DISTRIBUTION'
    cp.Incubation_Period_Gaussian_Mean = 3.0
    cp.Incubation_Period_Gaussian_Std_Dev = 1.0

    if(ID_STD > 0.01):
        cp.Infectious_Period_Distribution = 'GAMMA_DISTRIBUTION'
        cp.Infectious_Period_Scale = ID_MEAN*ID_STD*ID_STD
        cp.Infectious_Period_Shape = 1.0/ID_STD/ID_STD
    else:
        cp.Infectious_Period_Distribution = 'CONSTANT_DISTRIBUTION'
        cp.Infectious_Period_Constant = ID_MEAN

    cp.Enable_Infection_Rate_Overdispersion = 1
    cp.Enable_Infectivity_Scaling = 1

    cp.Enable_Disease_Mortality = 0

    cp.Acquisition_Transmission_Correlation = CORR_ACQ_TRANS

    # Network
    max_k = max_coeff_ref(NI_POWER)
    ni_coeff  = [np.exp(max_k[k1]+NI_LN_MULT[k1]) for k1 in range(len(max_k))]

    cp.Enable_Network_Infectivity = 1

    cp.Network_Infectivity_Coefficient = ni_coeff
    cp.Network_Infectivity_Exponent = NI_POWER
    cp.Network_Infectivity_Max_Export_Frac = NI_MAXFRAC
    cp.Network_Infectivity_Min_Connection = 1.0e-8
    cp.Network_Infectivity_Min_Distance = 1.0

    # Immunity
    cp.Enable_Immunity = 1
    cp.Enable_Immune_Decay = 0

    cp.Post_Infection_Acquisition_Multiplier = 0.0
    cp.Post_Infection_Transmission_Multiplier = 0.0
    cp.Post_Infection_Mortality_Multiplier = 0.0

    cp.Enable_Initial_Susceptibility_Distribution = 1
    cp.Susceptibility_Initialization_Distribution_Type = 'DISTRIBUTION_COMPLEX'

    cp.Maternal_Acquire_Config.Initial_Effect = 1.0
    cp.Maternal_Acquire_Config.Durability_Map.Times = [0, 50, 100, 150, 200, 250]
    cp.Maternal_Acquire_Config.Durability_Map.Values = [0.70, 0.57, 0.36, 0.17, 0.06, 0.00]

    # Interventions 
    cp.Enable_Interventions = 1
    cp.Campaign_Filename = gdata.camp_file

    # Adapted sampling
    cp.Individual_Sampling_Type = 'ADAPTED_SAMPLING_BY_IMMUNE_STATE'
    cp.Min_Node_Population_Samples =  100.0
    cp.Base_Individual_Sample_Rate = 1/AGENT_RATE
    cp.Relative_Sample_Rate_Immune = 0.05
    cp.Immune_Threshold_For_Downsampling = 1.0e-5
    cp.Immune_Downsample_Min_Age = 365.0

    # Multistrain

    num_strains = NOPV_BOXES + OPV_BOXES + 1
    log2_num_strains = np.ceil(np.log2(num_strains))

    cp.Enable_Strain_Tracking = 1
    cp.Enable_Genome_Dependent_Infectivity = 1
    cp.Enable_Genome_Mutation = 1
    cp.Enable_Label_By_Mutator = LABEL_MUTES

    cp.Number_of_Clades = 3
    cp.Log2_Number_of_Genomes_per_Clade = log2_num_strains

    list_multiply = np.ones(num_strains)
    list_multiply[:NOPV_BOXES] = NOPV_R0_MULT * R0_OPV / R0
    list_multiply[NOPV_BOXES:] = np.linspace(R0_OPV/R0, 1.0, num = OPV_BOXES + 1)
  
    cp.Genome_Infectivity_Multipliers = list_multiply.tolist()

    list_mutate = np.zeros(num_strains)
    list_mutate[:NOPV_BOXES] = NOPV_REV
    list_mutate[NOPV_BOXES:-1] = OPV_REV

    cp.Genome_Mutation_Rates = list_mutate.tolist()

    list_mlabel = np.zeros(num_strains)
    list_mlabel[-1] = 1

    config.parameters.Genome_Mutations_Labeled = list_mlabel.tolist()

    # Demographics
    cp.Enable_Demographics_Builtin = 0

    cp.Enable_Vital_Dynamics = 1

    cp.Enable_Birth = 1
    cp.Birth_Rate_Dependence = 'POPULATION_DEP_RATE'
    cp.Enable_Aging = 1
    cp.Age_Initialization_Distribution_Type = 'DISTRIBUTION_COMPLEX'
    cp.Enable_Natural_Mortality = 1
    enum_str = 'NONDISEASE_MORTALITY_BY_YEAR_AND_AGE_FOR_EACH_GENDER'
    cp.Death_Rate_Dependence = enum_str
    cp.Enable_Acquisition_Heterogeneity = 1

    cp.Demographics_Filenames = gdata.demog_files

    # HINT
    cp.Enable_Heterogeneous_Intranode_Transmission = 1

    # Reporting
    cp.Enable_Default_Reporting = 1
    cp.Enable_Demographics_Reporting = 1

    cp.Custom_Reports_Filename = REPORTS_FILE

    return config

# *****************************************************************************
