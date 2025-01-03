# *****************************************************************************
#
# *****************************************************************************

import json
import os

import global_data as gdata

import numpy as np

from emod_constants import O_FILE

# *****************************************************************************


def application(output_path):

    # Prep output dictionary
    SIM_IDX = gdata.sim_index
    key_str = '{:05d}'.format(SIM_IDX)
    parsed_dat = {key_str: dict()}

    # Calculate total attack rate and store in a json dict
    with open(os.path.join(output_path, 'InsetChart.json')) as fid01:
        inset_chart = json.load(fid01)

    ic_chan = inset_chart['Channels']
    new_inf = np.array(ic_chan['New Infections']['Data'])
    pop_vec = np.array(ic_chan['Statistical Population']['Data'])

    tot_pop = pop_vec[-1]
    max_inf = np.argmax(new_inf)
    tot_inf = np.sum(new_inf)
    epi_inf = np.sum(new_inf[:max_inf])

    parsed_dat[key_str]['atk_frac'] = int(tot_inf)/tot_pop
    parsed_dat[key_str]['herd_frac'] = int(epi_inf)/tot_pop

    # Write output dictionary
    with open(O_FILE, 'w') as fid01:
        json.dump(parsed_dat, fid01)

    return None

# *****************************************************************************
