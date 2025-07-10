# *****************************************************************************

import json
import os
import sys

import numpy as np
import scipy.optimize as spopt
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'local_python')))
from py_assets_common.emod_constants import EXP_V, NUM_SIMS, P_FILE, D_FILE

# *****************************************************************************

DIRNAME = 'experiment_covariance01'

# *****************************************************************************


def make_fig():

    # Figure setup
    fig01 = plt.figure(figsize=(8, 6))
    axs01 = fig01.add_subplot(111, label=None)
    plt.sca(axs01)

    axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
    axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
    axs01.set_axisbelow(True)

    axs01.set_xlabel(r'R$_{\mathrm{0}}$', fontsize=16)
    axs01.set_ylabel('Population Fraction', fontsize=16)

    axs01.set_xlim(0.50, 1.75)
    axs01.set_ylim(-0.01, 0.81)

    axs01.plot([1, 1], [0, 1], 'k--')
    axs01.plot([0, 2], [0, 0], 'k-')

    # Reference trajectory (Kermack-McKendric analytic solution)
    def KMlimt(x, R0):
        return 1-x-np.exp(-x*R0)

    xref = np.linspace(1.01, 2.0, 200)
    yref = np.zeros(xref.shape)
    for k1 in range(yref.shape[0]):
        yref[k1] = spopt.brentq(KMlimt, 1e-5, 1, args=(xref[k1]))
    axs01.plot(xref, yref, 'k-', lw=5.0,
               label='Attack Rate - Analytic Solution')

    # Sim outputs
    tpath = os.path.join('..', DIRNAME)

    with open(os.path.join(tpath, D_FILE)) as fid01:
        data_brick = json.load(fid01)

    with open(os.path.join(tpath, P_FILE)) as fid01:
        p_dict = json.load(fid01)

    N_SIMS = p_dict[NUM_SIMS]
    R0 = np.array(p_dict[EXP_V]['R0'])
    R0_VAR = np.round(np.array(p_dict[EXP_V]['R0_variance']), 2)
    ACQ_VAR = np.round(np.array(p_dict[EXP_V]['indiv_variance_acq']), 2)
    COR_VAL = np.round(np.array(p_dict[EXP_V]['correlation_acq_trans']), 2)
    INF_VAR = N_SIMS*[0.50]

    p_tuple = [(COR_VAL[k1], ACQ_VAR[k1], INF_VAR[k1], R0_VAR[k1])
               for k1 in range(N_SIMS)]
    p_tuple = sorted(list(set(p_tuple)))

    inf_mat = np.zeros(N_SIMS) - 1.0
    for sim_idx_str in data_brick:
        sim_idx = int(sim_idx_str)
        inf_mat[sim_idx] = np.array(data_brick[sim_idx_str]['atk_frac'])

    # Only successful sims will overwrite negatives with data
    fidx = (inf_mat >= 0.0)

    for tup in p_tuple:
        gidx = fidx & \
              (R0_VAR == tup[3]) & \
              (ACQ_VAR == tup[1]) & \
              (COR_VAL == tup[0])
        str_set = ['{:3.1f}'.format(tup[2]),
                   '{:3.1f}'.format(tup[1]),
                   '{:3.1f}'.format(tup[0])]
        label_str = 'Var$_{trans}$=' + str_set[0] + '; ' + \
                    'Var$_{acq}$=' + str_set[1] + '; ' + \
                    'Corr=' + str_set[2]
        axs01.plot(R0[gidx], inf_mat[gidx], lw=0.0,
                   marker='.', label=label_str)

    axs01.legend(fontsize=8)

    # Generate figures
    plt.tight_layout()
    plt.savefig(os.path.join('fig_attackrate01.png'))
    plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
