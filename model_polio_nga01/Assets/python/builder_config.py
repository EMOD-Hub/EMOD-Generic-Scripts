# *****************************************************************************
#
# Configuration file for simulation.
#
# *****************************************************************************

import global_data as gdata

import numpy as np

from emod_constants import CAMP_FILE, REPORTS_FILE, \
                           DIST_CMPLX, DIST_GAMMA, DIST_GAUSS

# *****************************************************************************


def max_coeff_ref(exp_vals):

    if (np.min(exp_vals) < 0.0 or np.max(exp_vals) > 8.0):
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
    R0_OPV_MULT = gdata.var_params['R0_OPV_mult']
    NOPV_R0_MULT = gdata.var_params['R0_nOPV_mult']
    NI_LN_MULT = gdata.var_params['net_inf_ln_mult']
    NI_POWER = gdata.var_params['net_inf_power']
    RUN_NUM = gdata.var_params['run_number']
    START_YEAR = gdata.var_params['start_year']
    RUN_YEARS = gdata.var_params['run_years']
    AGENT_RATE = gdata.var_params['agent_rate']

    # Config parameters object (read only dictionary)
    cp = config.parameters

    # Random number seed
    cp.Run_Number = RUN_NUM

    # Time
    cp.Start_Time = 365.0*(START_YEAR-gdata.base_year)
    cp.Simulation_Duration = 365.0*RUN_YEARS + gdata.t_step_days
    cp.Simulation_Timestep = gdata.t_step_days

    cp.Enable_Termination_On_Total_Wall_Time = 1
    cp.Wall_Time_Maximum_In_Minutes = gdata.max_clock

    # Intrahost
    ID_MEAN = gdata.inf_dur_mean
    ID_STD = gdata.inf_dur_std

    cp.Base_Infectivity_Distribution = DIST_GAMMA
    cp.Base_Infectivity_Shape = gdata.inf_tot_shape
    cp.Base_Infectivity_Scale = R0/ID_MEAN/gdata.inf_tot_shape

    cp.Incubation_Period_Distribution = DIST_GAUSS
    cp.Incubation_Period_Gaussian_Mean = 3.0
    cp.Incubation_Period_Gaussian_Std_Dev = 1.0

    cp.Infectious_Period_Distribution = DIST_GAMMA
    cp.Infectious_Period_Shape = ID_MEAN*ID_MEAN/ID_STD/ID_STD
    cp.Infectious_Period_Scale = ID_STD*ID_STD/ID_MEAN

    cp.Enable_Infection_Rate_Overdispersion = 1
    cp.Enable_Infectivity_Scaling = 1

    cp.Enable_Disease_Mortality = 0

    cp.Acquisition_Transmission_Correlation = 0.8

    # Network
    max_k = max_coeff_ref(NI_POWER)
    ni_coeff = [np.exp(max_k[k1]+NI_LN_MULT[k1]) for k1 in range(len(max_k))]

    cp.Enable_Network_Infectivity = 1

    cp.Network_Infectivity_Coefficient = ni_coeff
    cp.Network_Infectivity_Exponent = NI_POWER
    cp.Network_Infectivity_Max_Export_Frac = 0.05
    cp.Network_Infectivity_Min_Connection = 1.0e-7
    cp.Network_Infectivity_Min_Distance = 1.0

    # Immunity
    cp.Enable_Immunity = 1
    cp.Enable_Immune_Decay = 0

    cp.Post_Infection_Acquisition_Multiplier = 0.0
    cp.Post_Infection_Transmission_Multiplier = 0.0
    cp.Post_Infection_Mortality_Multiplier = 0.0

    cp.Enable_Initial_Susceptibility_Distribution = 1
    cp.Susceptibility_Initialization_Distribution_Type = DIST_CMPLX

    mac_x = [0, 50, 100, 150, 200, 250]
    mac_y = [0.70, 0.57, 0.36, 0.17, 0.06, 0.00]
    cp.Maternal_Acquire_Config.Initial_Effect = 1.0
    cp.Maternal_Acquire_Config.Durability_Map.Times = mac_x
    cp.Maternal_Acquire_Config.Durability_Map.Values = mac_y

    # Interventions
    cp.Enable_Interventions = 1
    cp.Campaign_Filename = CAMP_FILE

    # Adapted sampling
    cp.Individual_Sampling_Type = 'ADAPTED_SAMPLING_BY_IMMUNE_STATE'
    cp.Min_Node_Population_Samples = 100.0
    cp.Base_Individual_Sample_Rate = 1/AGENT_RATE
    cp.Relative_Sample_Rate_Immune = 0.05
    cp.Immune_Threshold_For_Downsampling = 1.0e-5
    cp.Immune_Downsample_Min_Age = 365.0

    # Multistrain
    CMOPV2 = gdata.boxes_sabin2
    CNOPV2 = gdata.boxes_nopv2
    num_strains = CNOPV2 + CMOPV2 + 1
    log2_num_strains = np.ceil(np.log2(num_strains))

    cp.Enable_Strain_Tracking = 1
    cp.Number_of_Clades = 2
    cp.Log2_Number_of_Genomes_per_Clade = log2_num_strains

    list_multiply = np.ones(num_strains)
    list_multiply[:CNOPV2] = NOPV_R0_MULT
    list_multiply[CNOPV2:] = np.linspace(R0_OPV_MULT, 1.0, num=CMOPV2+1)
    cp.Enable_Genome_Dependent_Infectivity = 1
    cp.Genome_Infectivity_Multipliers = list_multiply.tolist()

    list_mutate = np.zeros(num_strains)
    list_mutate[:CNOPV2] = gdata.rev_nopv2
    list_mutate[CNOPV2:-1] = gdata.rev_sabin2
    cp.Enable_Genome_Mutation = 1
    cp.Genome_Mutation_Rates = list_mutate.tolist()

    # Demographics
    cp.Enable_Demographics_Builtin = 0

    cp.Enable_Vital_Dynamics = 1

    cp.Enable_Birth = 1
    cp.Birth_Rate_Dependence = 'POPULATION_DEP_RATE'
    cp.Enable_Aging = 1
    cp.Age_Initialization_Distribution_Type = DIST_CMPLX
    cp.Enable_Natural_Mortality = 1
    enum_str = 'NONDISEASE_MORTALITY_BY_YEAR_AND_AGE_FOR_EACH_GENDER'
    cp.Death_Rate_Dependence = enum_str
    cp.Enable_Acquisition_Heterogeneity = 1

    cp.Demographics_Filenames = gdata.demog_files

    # Reporting
    cp.Enable_Default_Reporting = 1
    cp.Enable_Demographics_Reporting = 1

    cp.Custom_Reports_Filename = REPORTS_FILE

    # Logging
    cp.logLevel_StandardEventCoordinator = 'WARNING'

    # Memory
    cp.Memory_Usage_Halting_Threshold_Working_Set_MB = 15500
    cp.Memory_Usage_Warning_Threshold_Working_Set_MB = 15000

    return config

# *****************************************************************************
