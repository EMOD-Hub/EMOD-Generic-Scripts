# *****************************************************************************
#
# *****************************************************************************

import numpy as np

import matplotlib.patches as patch
import matplotlib.lines as line

from py_assets_common.emod_constants import POP_AGE_DAYS, CLR_M, CLR_F

# *****************************************************************************

PYR_TLOCS = [-12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12]
PYR_TLABS = [str(abs(val)) for val in PYR_TLOCS]

# *****************************************************************************


def haversine_dist(lat1, long1, lat2, long2):

    lat1 = np.pi*np.array(lat1)/180.0
    lat2 = np.pi*np.array(lat2)/180.0
    long1 = np.pi*np.array(long1)/180.0
    long2 = np.pi*np.array(long2)/180.0

    val1 = np.sin(0.5*lat2-0.5*lat1)**2
    val2 = np.cos(lat2)*np.cos(lat1)*np.sin(0.5*long2-0.5*long1)**2
    delt = 6371*2*np.arcsin(np.sqrt(val1+val2))

    return delt

# *****************************************************************************


def crs_proc(fert_file, XDAT, pyr_mat_avg, ss_demog=False):

    # Load UN WPP fertility distribution data
    fert_dat = np.loadtxt(fert_file, delimiter=',')
    fert_yr = fert_dat[0, :]
    fert_mat = fert_dat[1:, :]

    # Fertility distributions interpolated for years simulated
    fert_set = np.zeros((fert_mat.shape[0], XDAT.shape[0]))
    for k1 in range(fert_set.shape[0]):
        fert_set[k1, :] = np.interp(XDAT, fert_yr, fert_mat[k1, :])
    fert_set = fert_set/1000.0  # births/woman/year

    # Annual births implied by fertility distribution
    tot_fem = np.transpose((np.diff(pyr_mat_avg, axis=0)/2.0 +
                            pyr_mat_avg[:-1, :])/2.0)
    fertopt = fert_set
    if (ss_demog):
        fertopt = fert_set[:, 0]
        fertopt = fertopt[:, np.newaxis]
    fert_births = np.sum(fertopt*tot_fem, axis=0)

    # Calculate CRS probabilities
    crs_prob_vec = np.ones(fertopt.shape)       # P = 1
    crs_prob_vec = crs_prob_vec * 0.5           # P(female)
    crs_prob_vec = crs_prob_vec * fertopt       # P(gave birth during year)
    crs_prob_vec = crs_prob_vec * 9.0/12.0      # P(pregnant during year)
    crs_prob_vec = crs_prob_vec * (0.85*13/39 + 0.50*4/39 + 0.50*4/39)
                                                # P(infection leads to CRS)

    return (fert_births, crs_prob_vec)

# *****************************************************************************


def pyr_chart(axs01, pop_dat, pop_dat_err, yr_lab):

    axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
    axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
    axs01.set_axisbelow(True)

    axs01.set_xlabel('Percentage', fontsize=14)
    axs01.set_ylabel('Age (yrs)', fontsize=14)

    pdat_yr = np.array(POP_AGE_DAYS)/365.0

    axs01.set_xlim(min(PYR_TLOCS), max(PYR_TLOCS))
    axs01.set_ylim(pdat_yr[0], pdat_yr[-1])

    axs01.set_xticks(ticks=PYR_TLOCS)
    axs01.set_xticklabels(PYR_TLABS)

    ydat = pdat_yr - 2.5
    tpop = np.sum(pop_dat)

    pop_dat_n = 100*pop_dat/tpop
    pop_dat_n_err = 100*pop_dat_err/tpop

    axs01.barh(ydat[1:], pop_dat_n/2.0, height=4.75,
               xerr=pop_dat_n_err, color=CLR_F)
    axs01.barh(ydat[1:], -pop_dat_n/2.0, height=4.75,
               xerr=pop_dat_n_err, color=CLR_M)

    dval = 1e6
    tp_str = 'Total Pop\n{:5.1f}M'
    yr_str = '{:<4d}'
    if (tpop < 1e6):
        tp_str = 'Total Pop\n{:5.0f}k'
        dval = 1e3
    if (yr_lab < 1900):
        yr_str = 'Year {:<4d}'

    year_str = yr_str.format(yr_lab)
    axs01.text(-11, 92.5, year_str, fontsize=18)
    tpop_str = tp_str.format(tpop/dval)
    axs01.text(8, 87.5, tpop_str, fontsize=18, horizontalalignment='center')

    return None

# *****************************************************************************


def shape_patch(axs01, shp_pts, shp_brk, clr=[0.0, 0.5, 0.0], alph=1.0):

    skipnext = 0
    if (shp_brk[-1] != len(shp_pts)):
        shp_brk.append(len(shp_pts))

    for k1 in range(1, len(shp_brk)):
        if (skipnext > 0):
            skipnext = skipnext - 1
            continue

        subdat = np.array(shp_pts[shp_brk[k1-1]:shp_brk[k1]])
        cpos = k1

        while (cpos < len(shp_brk)-1):
            subdatadd = np.array(shp_pts[shp_brk[cpos]:shp_brk[cpos+1]])
            pn = np.sum((subdatadd[1:, 0]-subdatadd[:-1, 0]) *
                        (subdatadd[1:, 1]+subdatadd[:-1, 1]))
            if (pn < 0):
                subdat = np.vstack((subdat, subdatadd, subdat[-1, :]))
                skipnext = skipnext + 1
                cpos = cpos + 1
            else:
                break

        poly_shp = patch.Polygon(subdat, fc=clr, ec=clr,
                                 lw=1.0, ls='-', alpha=alph)
        axs01.add_patch(poly_shp)

    return None

# *****************************************************************************


def shape_line(axs01, shp_pts, shp_brk, clr=[0.0, 0.0, 0.0], wid=0.5):

    if (shp_brk[-1] != len(shp_pts)):
        shp_brk.append(len(shp_pts))

    for k1 in range(1, len(shp_brk)):
        subdat = np.array(shp_pts[shp_brk[k1-1]:shp_brk[k1]])
        line_shp = line.Line2D(subdat[:, 0], subdat[:, 1],
                               linewidth=wid, color=clr)
        axs01.add_line(line_shp)

    return None

# *****************************************************************************
