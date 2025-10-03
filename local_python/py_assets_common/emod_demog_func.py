# *****************************************************************************
#
# *****************************************************************************

import os

import numpy as np
import scipy.optimize as opt

from emod_api.demographics.demographics_overlay import DemographicsOverlay
from emod_api.demographics.age_distribution import AgeDistribution
from emod_api.demographics.mortality_distribution import MortalityDistribution
from emod_api.demographics.susceptibility_distribution import \
                                   SusceptibilityDistribution
from emod_api.demographics.PropertiesAndAttributes import \
                                   IndividualAttributes, NodeAttributes

from emod_api.demographics import DemographicsTemplates as DT

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
    age_y_list = (np.array(age_y_vec)/365.0).tolist()

    if (mort_age_in):
        mort_age_vec = mort_age_in
    else:
        mort_age_vec = MORT_XVAL
    mort_y_list = (np.array(mort_age_vec)/365.0).tolist()

    # Initial age and mortality distributions
    ind_age = AgeDistribution(ages_years=age_y_list,
                              cumulative_population_fraction=age_x)

    ind_mort = MortalityDistribution(ages_years=mort_y_list,
                                     calendar_years=mort_year,
                                     mortality_rate_matrix=mort_mat)

    ind_att = IndividualAttributes()
    ind_att.age_distribution = ind_age
    ind_att.mortality_distribution_male = ind_mort
    ind_att.mortality_distribution_female = ind_mort

    # Birth rate
    node_att = NodeAttributes()
    node_att.birth_rate = cb_rate

    # Overlay files
    dover_obj = DemographicsOverlay(idref=ref_name, nodes=node_list,
                                    individual_attributes=ind_att,
                                    node_attributes=node_att)

    nfname = DEMOG_FILE.rsplit('.', 1)[0] + '_vd{:04d}.json'.format(idx)
    nfname = os.path.join(PATH_OVERLAY, nfname)
    dover_obj.to_file(file_name=nfname)

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
    ind_sus = SusceptibilityDistribution()
    ind_sus.distribution_values = isus_x
    ind_sus.result_scale_factor = 1
    ind_sus.result_values = isus_y

    ind_att = IndividualAttributes()
    ind_att.susceptibility_distribution = ind_sus

    dover_obj = DemographicsOverlay(idref=ref_name, nodes=node_list,
                                    individual_attributes=ind_att)

    nfname = DEMOG_FILE.rsplit('.', 1)[0] + '_is{:04d}.json'.format(idx)
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

    nfname = DEMOG_FILE.rsplit('.', 1)[0] + '_is{:04d}.json'.format(idx)
    nfname = os.path.join(PATH_OVERLAY, nfname)
    dover_obj.to_file(file_name=nfname)

    return nfname

# *****************************************************************************


def demog_r0mult_over(ref_name, node_list, R0_mult, idx=0):

    if (not os.path.exists(PATH_OVERLAY)):
        os.mkdir(PATH_OVERLAY)

    dover_obj = DemographicsOverlay(idref=ref_name, nodes=node_list)

    nadict = dict()
    nadict['InfectivityMultiplier'] = R0_mult
    dover_obj.raw['Defaults']['NodeAttributes'].update(nadict)

    nfname = DEMOG_FILE.rsplit('.', 1)[0] + '_r0mult{:04d}.json'.format(idx)
    nfname = os.path.join(PATH_OVERLAY, nfname)
    dover_obj.to_file(file_name=nfname)

    return nfname

# *****************************************************************************


def demog_vd_calc(fname_pop, start_year, steady_state=False):

    # Load population data
    pop_input = np.loadtxt(fname_pop, dtype=int, delimiter=',')
    year_vec = pop_input[0, :]
    pop_mat = pop_input[1:, :] + 0.1

    # Estimate total initial population
    pop_init_vec = [np.interp(start_year, year_vec, pop_mat[idx, :])
                    for idx in range(pop_mat.shape[0])]

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
    brate_val = np.interp(start_year, year_vec[:-1], brate_vec)
    yrs_off = year_vec[:-1]-start_year
    yrs_dex = (yrs_off > 0)

    brmultx_01 = np.array([0.0] + (365.0*yrs_off[yrs_dex]).tolist())
    brmulty_01 = np.array([1.0] + (brate_vec[yrs_dex]/brate_val).tolist())
    brmultx_02 = np.zeros(2*len(brmultx_01)-1)
    brmulty_02 = np.zeros(2*len(brmulty_01)-1)

    brmultx_02[0::2] = brmultx_01[0:]
    brmulty_02[0::2] = brmulty_01[0:]
    brmultx_02[1::2] = brmultx_01[1:]-0.5
    brmulty_02[1::2] = brmulty_01[0:-1]

    if (steady_state):
        brmulty_02 = 1.0 + 0.0*brmulty_02

    brmx = brmultx_02.tolist()
    brmy = brmulty_02.tolist()

    pop_init = np.sum(pop_init_vec)
    age_init_cdf = np.cumsum(pop_init_vec[:-1])/pop_init
    age_x = [0] + age_init_cdf.tolist()

    b_rate = brate_val/365.0/1000.0
    mort_year = np.zeros(2*year_vec.shape[0]-3)

    mort_year[0::2] = year_vec[0:-1]
    mort_year[1::2] = year_vec[1:-1]-1e-3
    mort_year = mort_year.tolist()

    mort_mat = np.zeros((len(MORT_XVAL), len(mort_year)))

    mort_mat[0:-2:2, 0::2] = mortvecs
    mort_mat[1:-2:2, 0::2] = mortvecs
    mort_mat[0:-2:2, 1::2] = mortvecs[:, :-1]
    mort_mat[1:-2:2, 1::2] = mortvecs[:, :-1]
    mort_mat[-2:, :] = MAX_DAILY_MORT

    mort_mat_yr = 365.0*mort_mat
    mort_mat_yr = mort_mat_yr.tolist()

    if (steady_state):
        mort_vec = np.array([np.interp(start_year, mort_year, mort_mat[idx, :])
                             for idx in range(mort_mat.shape[0])])
        mort_year = [start_year]
        mort_mat_yr = 365.0*mort_vec[:, np.newaxis]
        mort_mat_yr = mort_mat_yr.tolist()
        mort_vec = mort_vec.tolist()

        forcing_vec = 12*[1.0]  # No seasonal forcing
        (_, age_x_eq, age_y_eq) = DT._computeAgeDist(b_rate, MORT_XVAL,
                                                     mort_vec, forcing_vec)
        age_x = (np.interp(POP_AGE_DAYS, age_y_eq, age_x_eq)).tolist()

    # Sanitize age cdf to insure it is increasing
    age_x[-1] = 1.0
    for k1 in range(2, len(age_x)):
        if (age_x[-k1] > (1.0 - k1*1e-6)):
            age_x[-k1] = (1.0 - k1*1e-6)
        else:
            break

    return (pop_init, mort_year, mort_mat_yr, age_x, b_rate, brmx, brmy)

# *****************************************************************************
