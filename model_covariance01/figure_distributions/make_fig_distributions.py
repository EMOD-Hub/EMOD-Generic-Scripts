# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'local_python')))
from py_assets_common.emod_constants import EXP_V, NUM_SIMS, P_FILE

# *****************************************************************************

DIRNAME = 'experiment_covariance01'

# *****************************************************************************


def make_fig():

    # Sim outputs
    tpath = os.path.join('..', DIRNAME)

    with open(os.path.join(tpath, P_FILE)) as fid01:
        p_dict = json.load(fid01)

    N_SIMS = p_dict[NUM_SIMS]
    R0_VAR = np.round(np.array(p_dict[EXP_V]['R0_variance']), 2)
    ACQ_VAR = np.round(np.array(p_dict[EXP_V]['indiv_variance_acq']), 2)
    COR_VAL = np.round(np.array(p_dict[EXP_V]['correlation_acq_trans']), 2)
    INF_VAR = N_SIMS*[0.50]

    p_tuple = [(COR_VAL[k1], ACQ_VAR[k1], INF_VAR[k1], R0_VAR[k1])
               for k1 in range(N_SIMS)]
    p_tuple = sorted(list(set(p_tuple)))

    # Figure
    fig01 = plt.figure(figsize=(12, 9))

    axs00 = fig01.add_subplot(1, 1, 1)
    axs00.spines['top'].set_color('none')
    axs00.spines['bottom'].set_color('none')
    axs00.spines['left'].set_color('none')
    axs00.spines['right'].set_color('none')
    axs00.tick_params(labelcolor='w', top=False, bottom=False,
                      left=False, right=False)

    for k1 in range(len(p_tuple)):
        axs01 = fig01.add_subplot(2, 2, k1+1, label=None)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)

        axs01.set_xlabel('Acquision Rate Multiplier')
        axs01.set_ylabel('Transmission Rate Multiplier')

        axs01.set_xlim(0, 6)
        axs01.set_ylim(0, 6)

        R0_VAR = p_tuple[k1][3]
        AQ_VAR = p_tuple[k1][1]
        rho = p_tuple[k1][0]

        R0_LN_SIG = np.sqrt(np.log(R0_VAR+1.0))
        R0_LN_MU = -0.5*R0_LN_SIG*R0_LN_SIG

        ACQ_LN_SIG = np.sqrt(np.log(AQ_VAR+1.0))
        ACQ_LN_MU = -0.5*ACQ_LN_SIG*ACQ_LN_SIG

        N = 1000

        risk_vec = np.random.lognormal(mean=ACQ_LN_MU,
                                       sigma=ACQ_LN_SIG, size=N)
        inf_vec = np.random.lognormal(mean=R0_LN_MU,
                                      sigma=R0_LN_SIG, size=N)
        corr_vec = inf_vec*(1 + rho*(risk_vec - 1))

        axs01.plot(risk_vec, corr_vec, marker='.',
                   lw=0.0, c='C{:d}'.format(k1))

    plt.tight_layout()
    plt.savefig('fig_distributions01.png')
    plt.close()

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
