# *****************************************************************************
#
# *****************************************************************************

import json
import os

import numpy as np
import scipy.optimize as opt

from emod_api.demographics.Demographics import DemographicsOverlay
from emod_api.demographics.PropertiesAndAttributes import IndividualAttributes

from emod_constants import DEMOG_FILE, PATH_OVERLAY, \
                           MORT_XVAL, POP_AGE_DAYS, MAX_DAILY_MORT

# *****************************************************************************


def demog_vd_over(ref_name, node_list, cb_rate,
                  mort_year, mort_mat, age_x, age_y=None, idx=0,
                  mort_age_in=None):

    if (not os.path.exists(PATH_OVERLAY)):
        os.mkdir(PATH_OVERLAY)

    if (age_y):
        age_y_vec = age_y
    else:
        age_y_vec = POP_AGE_DAYS

    if (mort_age_in):
        mort_age_vec = mort_age_in
    else:
        mort_age_vec = MORT_XVAL

    vd_over_dict = dict()

    vd_over_dict['Metadata'] = {'IdReference': ref_name}
    vd_over_dict['Defaults'] = {'IndividualAttributes': dict(),
                                'NodeAttributes': dict()}
    vd_over_dict['Nodes'] = [{'NodeID': nid.forced_id} for nid in node_list]

    vdodd = vd_over_dict['Defaults']
    vdodd['NodeAttributes'] = {'BirthRate': cb_rate}
    vdodd['IndividualAttributes'] = {'AgeDistribution': dict(),
                                     'MortalityDistributionMale': dict(),
                                     'MortalityDistributionFemale': dict()}

    vdoddiaad = vdodd['IndividualAttributes']['AgeDistribution']
    vdoddiaad['DistributionValues'] = [age_x]
    vdoddiaad['ResultScaleFactor'] = 1
    vdoddiaad['ResultValues'] = [age_y_vec]

    vdoddiamdm = vdodd['IndividualAttributes']['MortalityDistributionMale']
    vdoddiamdm['AxisNames'] = ['age', 'year']
    vdoddiamdm['AxisScaleFactors'] = [1, 1]
    vdoddiamdm['NumDistributionAxes'] = 2
    vdoddiamdm['NumPopulationGroups'] = [len(mort_age_vec), len(mort_year)]
    vdoddiamdm['PopulationGroups'] = [mort_age_vec, mort_year]
    vdoddiamdm['ResultScaleFactor'] = 1
    vdoddiamdm['ResultValues'] = mort_mat.tolist()

    vdoddiamdf = vdodd['IndividualAttributes']['MortalityDistributionFemale']
    vdoddiamdf['AxisNames'] = ['age', 'year']
    vdoddiamdf['AxisScaleFactors'] = [1, 1]
    vdoddiamdf['NumDistributionAxes'] = 2
    vdoddiamdf['NumPopulationGroups'] = [len(mort_age_vec), len(mort_year)]
    vdoddiamdf['PopulationGroups'] = [mort_age_vec, mort_year]
    vdoddiamdf['ResultScaleFactor'] = 1
    vdoddiamdf['ResultValues'] = mort_mat.tolist()

    nfname = DEMOG_FILE.rsplit('.', 1)[0] + '_vd{:03d}.json'.format(idx)
    nfname = os.path.join(PATH_OVERLAY, nfname)

    with open(nfname, 'w') as fid01:
        json.dump(vd_over_dict, fid01)

    return nfname

# *****************************************************************************


def min_fun(x1, age_year, age_prob, targ_frac):

    min_val = np.minimum(np.exp(x1*(age_year-0.65)), 1.0)*age_prob
    retval = np.sum(min_val)-targ_frac

    return retval

# *****************************************************************************


def demog_is_over(ref_name, node_list, R0, age_x, age_y=None, idx=0):

    if (not os.path.exists(PATH_OVERLAY)):
        os.mkdir(PATH_OVERLAY)

    if (age_y is None):
        age_y = POP_AGE_DAYS

    # Calculate initial susceptibilities
    targ_frac = 1.1*(1.0/R0)  # Tries to aim for Reff of 1.1

    # Implicit solve of exponential decay mapped onto age distribution. Target
    # area-under-the-curve is specified by targ_frac. Just aims to get close.
    # May break for very low target frac values (e.g., < 0.01)
    age_y_res = np.arange(1, 100*365, 30)
    age_x_res = np.interp(age_y_res, age_y, age_x)
    age_year = np.array(age_y_res[1:])/365.0
    age_prob = np.diff(np.array(age_x_res))

    arg_tup = (age_year, age_prob, targ_frac)
    iSP0 = opt.brentq(min_fun, a=-80, b=0, args=arg_tup)
    isus_x = [0] + (np.logspace(1.475, 4.540, 20, dtype=int)).tolist()
    isus_y = [round(np.minimum(np.exp(iSP0*(val/365.0-0.65)), 1.0), 4)
              for val in isus_x]

    # Initial susceptibility overlays
    ind_sus = IndividualAttributes.SusceptibilityDistribution()
    ind_sus.distribution_values = isus_x
    ind_sus.result_scale_factor = 1
    ind_sus.result_values = isus_y

    ind_att = IndividualAttributes()
    ind_att.susceptibility_distribution = ind_sus

    dover_obj = DemographicsOverlay(idref=ref_name, nodes=node_list,
                                    individual_attributes=ind_att)

    nfname = DEMOG_FILE.rsplit('.', 1)[0] + '_is{:03d}.json'.format(idx)
    nfname = os.path.join(PATH_OVERLAY, nfname)
    dover_obj.to_file(file_name=nfname)

    return nfname

# *****************************************************************************


def demog_is_over_precalc(ref_name, node_list, isus_x, isus_y, idx=0):

    if (not os.path.exists(PATH_OVERLAY)):
        os.mkdir(PATH_OVERLAY)

    # Initial susceptibility overlays
    ind_sus = IndividualAttributes.SusceptibilityDistribution()
    ind_sus.distribution_values = isus_x
    ind_sus.result_scale_factor = 1
    ind_sus.result_values = isus_y

    ind_att = IndividualAttributes()
    ind_att.susceptibility_distribution = ind_sus

    dover_obj = DemographicsOverlay(idref=ref_name, nodes=node_list,
                                    individual_attributes=ind_att)

    nfname = DEMOG_FILE.rsplit('.', 1)[0] + '_is{:03d}.json'.format(idx)
    nfname = os.path.join(PATH_OVERLAY, nfname)
    dover_obj.to_file(file_name=nfname)

    return nfname

# *****************************************************************************


def demog_ah_over(ref_name, node_list, R0_mult, ind_risk_var, idx=0):

    if (not os.path.exists(PATH_OVERLAY)):
        os.mkdir(PATH_OVERLAY)

    dover_obj = DemographicsOverlay(idref=ref_name, nodes=node_list)

    nadict = dict()
    nadict['InfectivityMultiplier'] = R0_mult
    dover_obj.raw['Defaults']['NodeAttributes'].update(nadict)

    iadict = dict()
    iadict['AcquisitionHeterogeneityVariance'] = ind_risk_var
    dover_obj.raw['Defaults']['IndividualAttributes'].update(iadict)

    nfname = DEMOG_FILE.rsplit('.', 1)[0] + '_ah{:03d}.json'.format(idx)
    nfname = os.path.join(PATH_OVERLAY, nfname)
    dover_obj.to_file(file_name=nfname)

    return nfname

# *****************************************************************************


def demog_vd_calc(year_vec, year_init, pop_mat, pop_init):

    # Calculate vital dynamics
    diff_ratio = (pop_mat[:-1, :-1]-pop_mat[1:, 1:])/pop_mat[:-1, :-1]
    t_delta = np.diff(year_vec)
    pow_vec = 365.0*t_delta
    mortvecs = 1.0-np.power(1.0-diff_ratio, 1.0/pow_vec)
    mortvecs = np.minimum(mortvecs, MAX_DAILY_MORT)
    mortvecs = np.maximum(mortvecs, 0.0)
    tot_pop = np.sum(pop_mat, axis=0)
    tpop_mid = (tot_pop[:-1]+tot_pop[1:])/2.0
    pop_corr = np.exp(-mortvecs[0, :]*pow_vec/2.0)

    brate_vec = np.round(pop_mat[0, 1:]/tpop_mid/t_delta*1000.0, 1)/pop_corr
    brate_val = np.interp(year_init, year_vec[:-1], brate_vec)
    yrs_off = year_vec[:-1]-year_init
    yrs_dex = (yrs_off > 0)

    brmultx_01 = np.array([0.0] + (365.0*yrs_off[yrs_dex]).tolist())
    brmulty_01 = np.array([1.0] + (brate_vec[yrs_dex]/brate_val).tolist())
    brmultx_02 = np.zeros(2*len(brmultx_01)-1)
    brmulty_02 = np.zeros(2*len(brmulty_01)-1)

    brmultx_02[0::2] = brmultx_01[0:]
    brmulty_02[0::2] = brmulty_01[0:]
    brmultx_02[1::2] = brmultx_01[1:]-0.5
    brmulty_02[1::2] = brmulty_01[0:-1]

    age_init_cdf = np.cumsum(pop_init[:-1])/np.sum(pop_init)
    age_x = [0] + age_init_cdf.tolist()

    birth_rate = brate_val/365.0/1000.0
    mort_year = np.zeros(2*year_vec.shape[0]-3)

    mort_year[0::2] = year_vec[0:-1]
    mort_year[1::2] = year_vec[1:-1]-1e-4
    mort_year = mort_year.tolist()

    mort_mat = np.zeros((len(MORT_XVAL), len(mort_year)))

    mort_mat[0:-2:2, 0::2] = mortvecs
    mort_mat[1:-2:2, 0::2] = mortvecs
    mort_mat[0:-2:2, 1::2] = mortvecs[:, :-1]
    mort_mat[1:-2:2, 1::2] = mortvecs[:, :-1]
    mort_mat[-2:, :] = MAX_DAILY_MORT

    return (mort_year, mort_mat, age_x, birth_rate, brmultx_02, brmulty_02)


# *****************************************************************************
