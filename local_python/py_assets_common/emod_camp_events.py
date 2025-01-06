# *****************************************************************************
#
# *****************************************************************************

from emod_api import schema_to_class as s2c
from emod_api.interventions import utils

from emod_constants import SPATH, YR_DAYS

# *****************************************************************************

CE = 'CampaignEvent'
SEC = 'StandardEventCoordinator'
CHWEC = 'CommunityHealthWorkerEventCoordinator'

NLHT = 'NodeLevelHealthTriggeredIV'
MID = 'MultiInterventionDistributor'
DI = 'DelayedIntervention'
VAX = 'Vaccine'
OI = 'OutbreakIndividual'

VEC_AGE = [0/12*365, 3/12*365, 5/12*365, 7/12*365, 9/12*365]
VEC_TAKE = [0.0, 0.0, 0.65, 0.92, 1.0]

# *****************************************************************************


def ce_import_pressure(node_list,
                       start_day=0.0, duration=1.0, magnitude=1.0,
                       age_yrs=40.0, clade=0, genome=0):

    # Import pressure
    camp_event = s2c.get_class_with_defaults(CE, SPATH)
    camp_coord = s2c.get_class_with_defaults(SEC, SPATH)
    camp_iv = s2c.get_class_with_defaults('ImportPressure', SPATH)

    node_set = utils.do_nodes(SPATH, node_list)

    camp_event.Event_Coordinator_Config = camp_coord
    camp_event.Start_Day = start_day
    camp_event.Nodeset_Config = node_set

    camp_coord.Intervention_Config = camp_iv

    camp_iv.Durations = [duration]
    camp_iv.Daily_Import_Pressures = [magnitude]
    camp_iv.Import_Age = age_yrs*YR_DAYS
    camp_iv.Import_Female_Prob = 0.0
    camp_iv.Import_Agent_MC_Weight = 1.0
    camp_iv.Clade = clade
    camp_iv.Genome = genome

    return camp_event

# *****************************************************************************


def ce_random_numbers(node_list, start_day=0.0, numbers=0):

    # Import pressure
    camp_event = s2c.get_class_with_defaults(CE, SPATH)
    camp_coord = s2c.get_class_with_defaults(SEC, SPATH)
    camp_iv = s2c.get_class_with_defaults('NodeRandomNumber', SPATH)

    node_set = utils.do_nodes(SPATH, node_list)

    camp_event.Event_Coordinator_Config = camp_coord
    camp_event.Start_Day = start_day
    camp_event.Nodeset_Config = node_set

    camp_coord.Intervention_Config = camp_iv

    camp_iv.Total_Numbers = numbers

    return camp_event

# *****************************************************************************


def ce_br_force(node_list, times, values,
                start_day=0.0):

    # Birth rate multiplier
    camp_event = s2c.get_class_with_defaults(CE, SPATH)
    camp_coord = s2c.get_class_with_defaults(SEC, SPATH)
    camp_iv = s2c.get_class_with_defaults('NodeBirthRateMult', SPATH)

    node_set = utils.do_nodes(SPATH, node_list)

    camp_event.Event_Coordinator_Config = camp_coord
    camp_event.Start_Day = start_day
    camp_event.Nodeset_Config = node_set

    camp_coord.Intervention_Config = camp_iv

    camp_iv.Multiplier_By_Duration.Times = times
    camp_iv.Multiplier_By_Duration.Values = values

    return camp_event

# *****************************************************************************


def ce_inf_force(node_list, step_init, step_width, step_size,
                 start_day=0.0):

    # Infectivity multiplier
    camp_event = s2c.get_class_with_defaults(CE, SPATH)
    camp_coord = s2c.get_class_with_defaults(SEC, SPATH)
    camp_iv = s2c.get_class_with_defaults('NodeInfectivityMult', SPATH)

    node_set = utils.do_nodes(SPATH, node_list)

    camp_event.Event_Coordinator_Config = camp_coord
    camp_event.Start_Day = start_day
    camp_event.Nodeset_Config = node_set

    camp_coord.Intervention_Config = camp_iv
    camp_coord.Number_Repetitions = -1  # Repeat forever
    camp_coord.Timesteps_Between_Repetitions = 365.0

    x_init = step_init % 365.0
    x_ends = (step_init+step_width) % 365.0
    y_bound = 1.0

    if (x_ends < x_init):
        y_bound = step_size

    xyvals = set([(0.0, y_bound),
                  (x_init-0.1, 1.0),
                  (x_init, step_size),
                  (x_ends-0.1, step_size),
                  (x_ends, 1.0),
                  (365.0, y_bound)])
    xyvals = sorted(list(xyvals), key=lambda val: val[0])

    camp_iv.Multiplier_By_Duration.Times = [val[0] for val in xyvals]
    camp_iv.Multiplier_By_Duration.Values = [val[1] for val in xyvals]

    return camp_event

# *****************************************************************************


def ce_RI(node_list,
          coverage=1.0, start_day=0.0, base_take=1.0,
          age_one=300.0, frac_two=None, age_two=475.0, age_std=90.0,
          acq_fact=0.0, age_dep=False):

    # Vaccine
    camp_event = s2c.get_class_with_defaults(CE, SPATH)
    camp_coord = s2c.get_class_with_defaults(SEC, SPATH)
    camp_iv01 = s2c.get_class_with_defaults(NLHT, SPATH)
    camp_iv02 = s2c.get_class_with_defaults(MID, SPATH)

    camp_iv03 = s2c.get_class_with_defaults(DI, SPATH)
    camp_iv04 = s2c.get_class_with_defaults(VAX, SPATH)

    node_set = utils.do_nodes(SPATH, node_list)

    camp_event.Event_Coordinator_Config = camp_coord
    camp_event.Start_Day = start_day
    camp_event.Nodeset_Config = node_set

    camp_coord.Intervention_Config = camp_iv01

    camp_iv01.Actual_IndividualIntervention_Config = camp_iv02
    camp_iv01.Demographic_Coverage = coverage
    camp_iv01.Trigger_Condition_List = ['Births']

    camp_iv02.Intervention_List = [camp_iv03]

    camp_iv03.Actual_IndividualIntervention_Configs = [camp_iv04]
    camp_iv03.Delay_Period_Distribution = "GAUSSIAN_DISTRIBUTION"
    camp_iv03.Delay_Period_Gaussian_Mean = age_one
    camp_iv03.Delay_Period_Gaussian_Std_Dev = age_std

    camp_iv04.Acquire_Config.Initial_Effect = 1.0
    camp_iv04.Vaccine_Take = base_take
    camp_iv04.Take_Reduced_By_Acquire_Immunity = acq_fact

    if (age_dep):
        camp_iv04.Take_By_Age_Multiplier.Times = VEC_AGE
        camp_iv04.Take_By_Age_Multiplier.Values = VEC_TAKE

    # Second RI dose
    if (frac_two):
        camp_iv05 = s2c.get_class_with_defaults(DI, SPATH)
        camp_iv06 = s2c.get_class_with_defaults(VAX, SPATH)

        camp_iv05.Actual_IndividualIntervention_Configs = [camp_iv06]
        camp_iv05.Coverage = frac_two
        camp_iv05.Delay_Period_Distribution = "GAUSSIAN_DISTRIBUTION"
        camp_iv05.Delay_Period_Gaussian_Mean = age_two
        camp_iv05.Delay_Period_Gaussian_Std_Dev = age_std

        camp_iv06.Acquire_Config.Initial_Effect = 1.0
        camp_iv06.Vaccine_Take = base_take
        camp_iv06.Take_Reduced_By_Acquire_Immunity = acq_fact

        camp_iv02.Intervention_List.append(camp_iv05)

        if (age_dep):
            camp_iv06.Take_By_Age_Multiplier.Times = VEC_AGE
            camp_iv06.Take_By_Age_Multiplier.Values = VEC_TAKE

    return camp_event

# *****************************************************************************


def ce_RI_OPV(node_list,
              coverage=1.0, start_day=0.0, base_take=0.7,
              age_one=180.0, age_std=15.0,
              clade=0, genome=0):

    # Outbreak Individual
    camp_event = s2c.get_class_with_defaults(CE, SPATH)
    camp_coord = s2c.get_class_with_defaults(SEC, SPATH)
    camp_iv01 = s2c.get_class_with_defaults(NLHT, SPATH)
    camp_iv02 = s2c.get_class_with_defaults(MID, SPATH)

    camp_iv03 = s2c.get_class_with_defaults(DI, SPATH)
    camp_iv04 = s2c.get_class_with_defaults(OI, SPATH)

    node_set = utils.do_nodes(SPATH, node_list)

    camp_event.Event_Coordinator_Config = camp_coord
    camp_event.Start_Day = start_day
    camp_event.Nodeset_Config = node_set

    camp_coord.Intervention_Config = camp_iv01

    camp_iv01.Actual_IndividualIntervention_Config = camp_iv02
    camp_iv01.Demographic_Coverage = coverage*base_take
    camp_iv01.Trigger_Condition_List = ['Births']

    camp_iv02.Intervention_List = [camp_iv03]

    camp_iv03.Actual_IndividualIntervention_Configs = [camp_iv04]
    camp_iv03.Delay_Period_Distribution = "GAUSSIAN_DISTRIBUTION"
    camp_iv03.Delay_Period_Gaussian_Mean = age_one
    camp_iv03.Delay_Period_Gaussian_Std_Dev = age_std

    camp_iv04.Clade = clade
    camp_iv04.Genome = genome
    camp_iv04.Ignore_Immunity = 0
    camp_iv04.Cost_To_Consumer = 1.0

    return camp_event

# *****************************************************************************


def ce_SIA(node_list, start_day=0.0, yrs_min=0.75, yrs_max=5.0,
           coverage=0.8, base_take=1.0, acq_fact=0.0, age_dep=False):

    # Vaccine
    camp_event = s2c.get_class_with_defaults(CE, SPATH)
    camp_coord = s2c.get_class_with_defaults(SEC, SPATH)
    camp_iv01 = s2c.get_class_with_defaults(DI, SPATH)
    camp_iv02 = s2c.get_class_with_defaults(VAX, SPATH)

    node_set = utils.do_nodes(SPATH, node_list)

    camp_event.Event_Coordinator_Config = camp_coord
    camp_event.Start_Day = start_day
    camp_event.Nodeset_Config = node_set

    camp_coord.Intervention_Config = camp_iv01
    camp_coord.Target_Demographic = 'ExplicitAgeRanges'
    camp_coord.Demographic_Coverage = coverage
    camp_coord.Target_Age_Min = yrs_min
    camp_coord.Target_Age_Max = yrs_max

    camp_iv01.Actual_IndividualIntervention_Configs = [camp_iv02]
    camp_iv01.Delay_Period_Distribution = "UNIFORM_DISTRIBUTION"
    camp_iv01.Delay_Period_Min = 0.0
    camp_iv01.Delay_Period_Max = 14.0

    camp_iv02.Acquire_Config.Initial_Effect = 1.0
    camp_iv02.Vaccine_Take = base_take
    camp_iv02.Take_Reduced_By_Acquire_Immunity = acq_fact

    if (age_dep):
        camp_iv02.Take_By_Age_Multiplier.Times = VEC_AGE
        camp_iv02.Take_By_Age_Multiplier.Values = VEC_TAKE

    return camp_event

# *****************************************************************************


def ce_vax_AMT(node_list,
               start_day=0.0, only_group=None, coverage=1.0,
               acq_eff=0.0, mrt_eff=0.0, trn_eff=0.0):

    # Vaccine
    camp_event = s2c.get_class_with_defaults(CE, SPATH)
    camp_coord = s2c.get_class_with_defaults(SEC, SPATH)
    camp_iv = s2c.get_class_with_defaults(VAX, SPATH)

    node_set = utils.do_nodes(SPATH, node_list)

    camp_event.Event_Coordinator_Config = camp_coord
    camp_event.Start_Day = start_day
    camp_event.Nodeset_Config = node_set

    camp_coord.Intervention_Config = camp_iv

    if (only_group):
        camp_coord.Property_Restrictions_Within_Node = only_group

    camp_iv.Vaccine_Take = coverage
    camp_iv.Acquire_Config.Initial_Effect = acq_eff
    camp_iv.Mortality_Config.Initial_Effect = mrt_eff
    camp_iv.Transmit_Config.Initial_Effect = trn_eff

    return camp_event

# *****************************************************************************


def ce_OPV_SIA(node_list,
               start_day=0.0, coverage=1.0, take=0.70, clade=0, genome=0):

    # OutbreakIndividual
    camp_event = s2c.get_class_with_defaults(CE, SPATH)
    camp_coord = s2c.get_class_with_defaults(SEC, SPATH)
    camp_iv = s2c.get_class_with_defaults(OI, SPATH)

    node_set = utils.do_nodes(SPATH, node_list)

    camp_event.Event_Coordinator_Config = camp_coord
    camp_event.Start_Day = start_day
    camp_event.Nodeset_Config = node_set

    camp_coord.Intervention_Config = camp_iv
    camp_coord.Demographic_Coverage = coverage*take
    camp_coord.Target_Demographic = 'ExplicitAgeRanges'
    camp_coord.Target_Age_Min = 0.50
    camp_coord.Target_Age_Max = 5.00

    camp_iv.Clade = clade
    camp_iv.Genome = genome
    camp_iv.Ignore_Immunity = 0
    camp_iv.Cost_To_Consumer = 1.0

    return camp_event

# *****************************************************************************


def ce_quarantine(node_list, trigger,
                  start_day=0.0, coverage=1.0, delay=0.0, effect=0.0):

    # Vaccine
    camp_event = s2c.get_class_with_defaults(CE, SPATH)
    camp_coord = s2c.get_class_with_defaults(CHWEC, SPATH)
    camp_iv = s2c.get_class_with_defaults(VAX, SPATH)

    node_set = utils.do_nodes(SPATH, node_list)

    camp_event.Event_Coordinator_Config = camp_coord
    camp_event.Start_Day = start_day
    camp_event.Nodeset_Config = node_set

    camp_coord.Intervention_Config = camp_iv
    camp_coord.Trigger_Condition_List = [trigger]
    camp_coord.Duration = 1000
    camp_coord.Waiting_Period = 1000
    camp_coord.Max_Distributed_Per_Day = 1e9
    camp_coord.Days_Between_Shipments = 1000
    camp_coord.Amount_In_Shipment = 0
    camp_coord.Initial_Amount_Distribution = 'CONSTANT_DISTRIBUTION'
    camp_coord.Initial_Amount_Constant = 1e9

    camp_iv.Vaccine_Take = coverage

    camp_iv.Transmit_Config.Initial_Effect = 1.0
    camp_iv.Transmit_Config.Durability_Map.Times = [0.0, delay]
    camp_iv.Transmit_Config.Durability_Map.Values = [0.0, effect]

    return camp_event

# *****************************************************************************


def ce_matrix_swap(node_list, prop_name, matrix,
                   start_day=0.0):

    # ChangeIPMatrix
    camp_event = s2c.get_class_with_defaults(CE, SPATH)
    camp_coord = s2c.get_class_with_defaults(SEC, SPATH)
    camp_iv = s2c.get_class_with_defaults('ChangeIPMatrix', SPATH)

    node_set = utils.do_nodes(SPATH, node_list)

    camp_event.Event_Coordinator_Config = camp_coord
    camp_event.Start_Day = start_day
    camp_event.Nodeset_Config = node_set

    camp_coord.Intervention_Config = camp_iv

    camp_iv.Property_Name = prop_name
    camp_iv.New_Matrix = matrix

    return camp_event

# *****************************************************************************


def ce_visit_nodes(node_orig, node_dest,
                   start_day=0.0, only_group=None, fraction=0.0, duration=1.0):

    # MigrateIndividuals
    camp_event = s2c.get_class_with_defaults(CE, SPATH)
    camp_coord = s2c.get_class_with_defaults(SEC, SPATH)
    camp_iv = s2c.get_class_with_defaults('MigrateIndividuals', SPATH)

    node_set = utils.do_nodes(SPATH, [node_orig])

    camp_event.Event_Coordinator_Config = camp_coord
    camp_event.Start_Day = start_day
    camp_event.Nodeset_Config = node_set

    camp_coord.Intervention_Config = camp_iv
    camp_coord.Demographic_Coverage = fraction

    if (only_group):
        camp_coord.Property_Restrictions_Within_Node = only_group

    camp_iv.NodeID_To_Migrate_To = node_orig
    camp_iv.Duration_Before_Leaving_Distribution = 'CONSTANT_DISTRIBUTION'
    camp_iv.Duration_Before_Leaving_Constant = 0
    camp_iv.Duration_At_Node_Distribution = 'CONSTANT_DISTRIBUTION'
    camp_iv.Duration_At_Node_Constant = duration

    return camp_event

# *****************************************************************************
