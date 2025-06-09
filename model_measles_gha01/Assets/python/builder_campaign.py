#********************************************************************************
#
# Builds a campaign file for input to the DTK.
#
#********************************************************************************

import json
import os

import global_data as gdata

import numpy as np

import emod_api.campaign as camp_module

from emod_camp_events import ce_inf_force

from emod_constants import CAMP_FILE

#********************************************************************************


def campaignBuilder():

    # Variables for this simulation
    TIME_START   = gdata.start_time
    PEAK_SIZE    = gdata.var_params['R0_peak_magnitude']
    PEAK_TIME    = gdata.var_params['R0_peak_day']
    PEAK_WIDE    = gdata.var_params['R0_peak_width']

    # Note: campaign module itself is the file object; no Campaign class
    ALL_NODES = gdata.demog_object.node_ids
    CAMP_FILENAME =  'campaign.json'






    # ***** Events *****
    node_dict = gdata.demog_node
    node_opts = list(node_dict.keys())

    # Add MCV1 RI
    mcv1_dat = np.loadtxt(os.path.join('Assets','data', 'GHA_MCV1.csv'), delimiter=',')
    with open(os.path.join('Assets','data', 'GHA_MCV1.json')) as fid01:
        mcv1_dict = json.load(fid01)

    time_vec  = np.array(mcv1_dict['timevec']) - TIME_START
    for k1 in range(len(mcv1_dict['namevec'])):
        reg_name  = mcv1_dict['namevec'][k1]
        mcv1_vec  = mcv1_dat[k1,:]
        node_list = list()

        for node_name in node_opts:
            if((node_name == reg_name) or (node_name.startswith(reg_name+':'))):
                node_list.append(node_dict[node_name])

        if(not node_list):
            continue

        if(np.amin(time_vec) <= 0.0):
            init_mcv1 = np.interp(0.0, time_vec, mcv1_vec)
        else:
            init_mcv1 = np.mean(mcv1_vec[:3])

        time_list = [0.0]       + (time_vec[time_vec>0.0]).tolist() + [365.0*100]
        mcv1_list = [init_mcv1] + (mcv1_vec[time_vec>0.0]).tolist() + [np.mean(mcv1_vec[-3:])]

        pdict = {'startday':       TIME_START ,
                 'nodes':          node_list  ,
                 'x_vals':         time_list   ,
                 'y_vals':         mcv1_list   }

        camp_module.add(IV_MCV1(pdict))


    # Add MCV SIAs
    with open(os.path.join('Assets','data','GHA_MCV_SIA.json')) as fid01:
        dict_sia = json.load(fid01)

    for sia_name in dict_sia:
        sia_obj = dict_sia[sia_name]
        start_val  = sia_obj['date']
        if(start_val < TIME_START):
            continue

        SIA_COVER  = gdata.var_params['SIA_cover_{:s}'.format(sia_name)]
        age_min    = sia_obj['age_min']
        age_max    = sia_obj['age_max']
        targ_frac  = sia_obj['targ_frac']
        node_list  = list()
        for targ_val in sia_obj['nodes']:
            for node_name in node_opts:
                if((node_name == targ_val) or (node_name.startswith(targ_val+':'))):
                    node_list.append(node_dict[node_name])

        if(not node_list):
           continue

        pdict = {'startday': start_val ,
                 'nodes': node_list ,
                 'agemin': age_min ,
                 'agemax': age_max ,
                 'coverage': SIA_COVER*targ_frac }

        camp_module.add(IV_SIA(pdict))







    # Add seasonality
    start_day = 365.0*(gdata.start_year-gdata.base_year)
    camp_event = ce_inf_force(ALL_NODES, PEAK_TIME, PEAK_WIDE, PEAK_SIZE,
                              start_day=start_day, dt=gdata.t_step_days)
    camp_module.add(camp_event)

    # Add infectivity trough
    start_day = 365.0*(2020.0-gdata.base_year)
    camp_event = ce_inf_force(ALL_NODES, 0.0, 365.0*2.0, 0.75, nreps=1
                              start_day=start_day, dt=gdata.t_step_days)
    camp_module.add(camp_event)

    # End file construction
    camp_module.save(filename=CAMP_FILE)

    return None

# *****************************************************************************














# Routine immunization for MCV1
def IV_MCV1(params=dict()):

  SCHEMA_PATH   =  gdata.schema_path

  camp_event = s2c.get_class_with_defaults('CampaignEvent',                            SCHEMA_PATH)
  camp_coord = s2c.get_class_with_defaults('StandardEventCoordinator',                 SCHEMA_PATH)
  camp_iv01  = s2c.get_class_with_defaults('NodeLevelHealthTriggeredIVScaleUpSwitch',  SCHEMA_PATH)
  camp_iv02  = s2c.get_class_with_defaults('DelayedIntervention',                      SCHEMA_PATH)
  camp_iv03  = s2c.get_class_with_defaults('Vaccine',                                  SCHEMA_PATH)
  camp_wane  = s2c.get_class_with_defaults('WaningEffect',                             SCHEMA_PATH)

  node_set   = utils.do_nodes(SCHEMA_PATH, params['nodes'])

  camp_event.Event_Coordinator_Config                   = camp_coord
  camp_event.Start_Day                                  = params['startday']
  camp_event.Nodeset_Config                             = node_set

  camp_coord.Intervention_Config                        = camp_iv01

  camp_iv01.Actual_IndividualIntervention_Config        = camp_iv02
  camp_iv01.Demographic_Coverage                        = 1.0                 # Required, not used
  camp_iv01.Trigger_Condition_List                      = ['Births']
  camp_iv01.Demographic_Coverage_Time_Profile           = 'InterpolationMap'
  camp_iv01.Coverage_vs_Time_Interpolation_Map.Times    = params['x_vals']
  camp_iv01.Coverage_vs_Time_Interpolation_Map.Values   = params['y_vals']
  camp_iv01.Not_Covered_IndividualIntervention_Configs  = []                  # Breaks if not present

  camp_iv02.Actual_IndividualIntervention_Configs       = [camp_iv03]
  camp_iv02.Delay_Period_Distribution                   = "GAUSSIAN_DISTRIBUTION"
  camp_iv02.Delay_Period_Gaussian_Mean                  =  300.0
  camp_iv02.Delay_Period_Gaussian_Std_Dev               =   90.0

  camp_iv03.Acquire_Config                              = camp_wane

  camp_wane.Initial_Effect                              =    1.0

  return camp_event

#********************************************************************************

# SIAs for MCV
def IV_SIA(params=dict()):

  SCHEMA_PATH   =  gdata.schema_path

  camp_event = s2c.get_class_with_defaults('CampaignEvent',             SCHEMA_PATH)
  camp_coord = s2c.get_class_with_defaults('StandardEventCoordinator',  SCHEMA_PATH)
  camp_iv01  = s2c.get_class_with_defaults('Vaccine',                   SCHEMA_PATH)
  camp_wane  = s2c.get_class_with_defaults('WaningEffect',              SCHEMA_PATH)

  node_set   = utils.do_nodes(SCHEMA_PATH, params['nodes'])

  camp_event.Event_Coordinator_Config       = camp_coord
  camp_event.Start_Day                      = params['startday']
  camp_event.Nodeset_Config                 = node_set

  camp_coord.Intervention_Config            = camp_iv01
  camp_coord.Target_Demographic             = 'ExplicitAgeRanges'
  camp_coord.Demographic_Coverage           = params['coverage']
  camp_coord.Target_Age_Min                 = params['agemin']/365.0
  camp_coord.Target_Age_Max                 = params['agemax']/365.0

  camp_iv01.Acquire_Config                  = camp_wane

  camp_wane.Initial_Effect                  = 1.0

  return camp_event

#********************************************************************************
