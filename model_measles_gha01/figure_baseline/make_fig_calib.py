# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))

from py_assets_common.emod_constants import NUM_SIMS, P_FILE, D_FILE, \
                                            EXP_C, EXP_V

# *****************************************************************************

DIRNAMES = ['experiment_gha_base01']

# *****************************************************************************


def make_fig():

    for dirname in DIRNAMES:

        # Sim outputs
        tpath = os.path.join('..', dirname)

        with open(os.path.join(tpath, P_FILE)) as fid01:
            param_dict = json.load(fid01)

        with open(os.path.join(tpath, D_FILE)) as fid01:
            dbrick = json.load(fid01)

        nsims = int(param_dict[NUM_SIMS])
        tvals = dbrick.pop('tstamps')

        scale_vec = np.zeros((nsims, 1))-1
        cal_vec = np.zeros(nsims)

        for sim_idx_str in dbrick:
            if (not sim_idx_str.isdigit()):
                continue
            sim_idx = int(sim_idx_str)
            sim_obj = dbrick[sim_idx_str]

            scale_vec[sim_idx, 0] = sim_obj['rep_rate']
            cal_vec[sim_idx] = sim_obj['cal_val']

        gidx = (scale_vec[:,0] >= 0)

        sparam = param_dict[EXP_V]
        skeys = list(sparam.keys())
        skeys.remove('run_number')
        print(skeys)
        nfigs = len(skeys)

        # Figure
        fig01 = plt.figure(figsize=(8*nfigs, 8*nfigs))

        for k1 in range(nfigs):
            for k2 in range(k1, nfigs):

                if (k1==k2):
                    continue

                axs01 = fig01.add_subplot(nfigs, nfigs, k1*nfigs+k2+1)
                axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
                axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
                axs01.set_axisbelow(True)

                xval = np.array(sparam[skeys[k1]])
                yval = np.array(sparam[skeys[k2]])

                axs01.set_xlabel(skeys[k1])
                axs01.set_ylabel(skeys[k2])
                axs01.scatter(xval[gidx], yval[gidx], c=cal_vec[gidx],
                              vmin=-4.0e3, vmax=-2.5e3)


        print(min(cal_vec[gidx]), max(cal_vec[gidx]))
        plt.tight_layout()
        plt.savefig('fig_baseline01.png')
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
