# *****************************************************************************
#
# Configuration file for simulation.
#
# *****************************************************************************

import global_data as gdata

from emod_constants import CAMP_FILE, REPORTS_FILE

# *****************************************************************************


def update_config_obj(config):

    # Variables for this simulation
    RUN_NUM = gdata.var_params['run_number']
    R0 = gdata.var_params['R0']
    MAT_FACTOR = gdata.var_params['mat_factor']

    # Config parameters object (read only dictionary)
    cp = config.parameters

    # Random number seed
    cp.Run_Number = RUN_NUM

    # Time
    cp.Start_Time = 365.0*(gdata.start_year-gdata.base_year)
    cp.Simulation_Duration = 365.0*gdata.run_years
    cp.Simulation_Timestep = gdata.t_step_days

    cp.Enable_Termination_On_Total_Wall_Time = 1
    cp.Wall_Time_Maximum_In_Minutes = gdata.max_clock

    # Intrahost
    cp.Base_Infectivity_Distribution = 'GAMMA_DISTRIBUTION'
    cp.Base_Infectivity_Scale = R0/gdata.inf_prd_mean
    cp.Base_Infectivity_Shape = 1.0

    cp.Incubation_Period_Distribution = 'GAMMA_DISTRIBUTION'
    cp.Incubation_Period_Scale = 1.0
    cp.Incubation_Period_Shape = 3.5

    cp.Infectious_Period_Distribution = 'GAUSSIAN_DISTRIBUTION'
    cp.Infectious_Period_Gaussian_Mean = gdata.inf_prd_mean
    cp.Infectious_Period_Gaussian_Std_Dev = 2.0

    cp.Enable_Nonuniform_Shedding = 1.0
    cp.Shedding_Distribution_Alpha = 10.0
    cp.Shedding_Distribution_Beta = 10.0

    cp.Enable_Disease_Mortality = 0

    cp.Symptomatic_Infectious_Offset = 11.0

    # Immunity
    cp.Enable_Immunity = 1
    cp.Enable_Immune_Decay = 0

    cp.Post_Infection_Acquisition_Multiplier = 0.0
    cp.Post_Infection_Mortality_Multiplier = 0.0
    cp.Post_Infection_Transmission_Multiplier = 0.0

    cp.Maternal_Acquire_Config.Initial_Effect = MAT_FACTOR
    cp.Maternal_Acquire_Config.Enable_Box_Duration_Distribution = 1
    enum_str = 'GAUSSIAN_DISTRIBUTION'
    cp.Maternal_Acquire_Config.Box_Duration_Distribution = enum_str
    cp.Maternal_Acquire_Config.Box_Duration_Gaussian_Mean = 120.0
    cp.Maternal_Acquire_Config.Box_Duration_Gaussian_Std_Dev = 67.0

    cp.Enable_Initial_Susceptibility_Distribution = 1
    cp.Susceptibility_Initialization_Distribution_Type = 'DISTRIBUTION_COMPLEX'

    # Interventions
    cp.Enable_Interventions = 1
    cp.Campaign_Filename = CAMP_FILE

    # Infectivity
    cp.Enable_Acquisition_Heterogeneity = 0
    cp.Enable_Infection_Rate_Overdispersion = 0
    cp.Enable_Infectivity_Reservoir = 1

    # Adapted sampling
    cp.Individual_Sampling_Type = 'ADAPTED_SAMPLING_BY_IMMUNE_STATE'
    cp.Min_Node_Population_Samples = 20.0
    cp.Base_Individual_Sample_Rate = 1.0
    cp.Relative_Sample_Rate_Immune = 0.01
    cp.Immune_Threshold_For_Downsampling = 1.0e-5
    cp.Immune_Downsample_Min_Age = 365.0

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

    cp.Demographics_Filenames = gdata.demog_files

    # Reporting
    cp.Enable_Default_Reporting = 1
    cp.Enable_Demographics_Reporting = 1
    cp.Enable_Event_DB = 1
    cp.SQL_Events = ["NewInfection"]

    cp.Custom_Reports_Filename = REPORTS_FILE

    return config

# *****************************************************************************
