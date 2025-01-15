# *****************************************************************************
#
# *****************************************************************************

import os

import global_data as gdata

import numpy as np

from emod_constants import RST_FILE, RST_TIME, RST_NODE, RST_CLADE, \
                           RST_GENOME, RST_TOT_INF, RST_NEW_INF


# *****************************************************************************


def application(timestep):

    time_val = float(timestep)

    # Check for early sim abort
    if (gdata.inproc_abort_bool and time_val > gdata.inproc_abort_time):

        # Only check once
        gdata.inproc_abort_bool = False
        print('Checking contagion at time {:.0f}'.format(time_val))

        infdat = np.loadtxt(os.path.join('output', RST_FILE),
                            delimiter=',', skiprows=1, ndmin=2)

        sabin_clade = 0
        cVDPV_genome = gdata.boxes_sabin2 + gdata.boxes_nopv2
        targ_time = infdat[-1, RST_TIME]
        print('Last timestep in report: {:.0f}'.format(targ_time))

        # Sabin derived cVDPV2 infections in most recent timestep
        gidx = (infdat[:, RST_CLADE] == sabin_clade)
        gidx = (infdat[:, RST_GENOME] == cVDPV_genome) & gidx
        gidx = (infdat[:, RST_TIME] == targ_time) & gidx

        inf_sum = np.sum(infdat[gidx, RST_TOT_INF] + infdat[gidx, RST_NEW_INF])
        if (inf_sum == 0):
            return "ABORT"

    return None

# *****************************************************************************
