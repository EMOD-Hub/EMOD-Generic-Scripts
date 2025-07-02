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
    NI_LN_MULT = gdata.var_params['net_inf_ln_mult']
    NI_POWER = gdata.var_params['net_inf_power']
    NI_MAXFRAC = gdata.var_params['net_inf_maxfrac']
    RUN_NUM = gdata.var_params['run_number']
    END_YEAR = gdata.var_params['end_year']
    AGENT_RATE = gdata.var_params['agent_rate']
    CORR_ACQ_TRANS = gdata.var_params['corr_acq_trans']

    # Config parameters object (read only dictionary)
    cp = config.parameters

    # Random number seed
    cp.Run_Number = RUN_NUM

    # Time
    RUN_YEARS = END_YEAR - gdata.start_year
    cp.Start_Time = 365.0*(gdata.start_year-gdata.base_year)
    cp.Simulation_Duration = 365.0*RUN_YEARS + gdata.t_step_days
    cp.Simulation_Timestep = gdata.t_step_days

    cp.Enable_Termination_On_Total_Wall_Time = 1
    cp.Wall_Time_Maximum_In_Minutes = gdata.max_clock

    # Intrahost
    cp.Base_Infectivity_Distribution = DIST_GAMMA
    cp.Base_Infectivity_Scale = R0/gdata.inf_dur_mean
    cp.Base_Infectivity_Shape = 1.0

    cp.Incubation_Period_Distribution = DIST_GAMMA
    cp.Incubation_Period_Scale = 1.0
    cp.Incubation_Period_Shape = 3.5

    cp.Infectious_Period_Distribution = DIST_GAUSS
    cp.Infectious_Period_Gaussian_Mean = gdata.inf_dur_mean
    cp.Infectious_Period_Gaussian_Std_Dev = 2.0

    cp.Enable_Nonuniform_Shedding = 1.0
    cp.Shedding_Distribution_Alpha = 10.0
    cp.Shedding_Distribution_Beta = 10.0

    cp.Enable_Infection_Rate_Overdispersion = 1
    cp.Enable_Infectivity_Scaling = 1

    cp.Symptomatic_Infectious_Offset = 11.0

    cp.Enable_Disease_Mortality = 0

    cp.Acquisition_Transmission_Correlation = CORR_ACQ_TRANS

    # Immunity
    cp.Enable_Immunity = 1
    cp.Enable_Immune_Decay = 0

    cp.Post_Infection_Acquisition_Multiplier = 0.0
    cp.Post_Infection_Transmission_Multiplier = 0.0
    cp.Post_Infection_Mortality_Multiplier = 0.0

    cp.Maternal_Acquire_Config.Initial_Effect = 1.0
    cp.Maternal_Acquire_Config.Enable_Box_Duration_Distribution = 1
    cp.Maternal_Acquire_Config.Box_Duration_Distribution = DIST_GAUSS
    cp.Maternal_Acquire_Config.Box_Duration_Gaussian_Mean = 150.0
    cp.Maternal_Acquire_Config.Box_Duration_Gaussian_Std_Dev = 80.0

    cp.Enable_Initial_Susceptibility_Distribution = 1
    cp.Susceptibility_Initialization_Distribution_Type = DIST_CMPLX

    # Interventions
    cp.Enable_Interventions = 1
    cp.Campaign_Filename = CAMP_FILE

    # Network Infectivity
    ni_coeff = np.exp(max_coeff_ref(NI_POWER) + NI_LN_MULT)

    cp.Enable_Network_Infectivity = 1

    cp.Network_Infectivity_Coefficient = [ni_coeff]
    cp.Network_Infectivity_Exponent = [NI_POWER]
    cp.Network_Infectivity_Max_Export_Frac = NI_MAXFRAC
    cp.Network_Infectivity_Min_Distance = 1

    cp.Enable_Infectivity_Reservoir = 1

    # Adapted sampling
    cp.Individual_Sampling_Type = 'ADAPTED_SAMPLING_BY_IMMUNE_STATE'
    cp.Min_Node_Population_Samples = gdata.demog_min_pop
    cp.Base_Individual_Sample_Rate = 1.0/AGENT_RATE
    cp.Relative_Sample_Rate_Immune = 0.02
    cp.Immune_Threshold_For_Downsampling = 1.0e-5
    cp.Immune_Downsample_Min_Age = 365.0

    # Demographic parameters
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
    cp.Enable_Event_DB = 1
    cp.SQL_Start_Time = 365.0*(gdata.start_year_log-gdata.base_year)
    cp.SQL_Events = ["NewlySymptomatic"]

    cp.Enable_Spatial_Output = 0

    cp.Custom_Reports_Filename = REPORTS_FILE

    # Logging
    cp.logLevel_StandardEventCoordinator = 'WARNING'

    # Memory
    cp.Memory_Usage_Halting_Threshold_Working_Set_MB = 15500
    cp.Memory_Usage_Warning_Threshold_Working_Set_MB = 15000

    return config

# *****************************************************************************
