# *****************************************************************************

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join('..', 'refdat_namesets')))
sys.path.insert(0, os.path.abspath(os.path.join('..', 'local_python')))
from aux_namematch import reprule

# *****************************************************************************


def make_pop_dat(TLC=''):

    # Name references
    tlc_wpp_dict = dict()

    fname = os.path.join('..', 'refdat_namesets', 'tlc_wpp_countries.json')
    with open(fname) as fid01:
        n_dict = json.load(fid01)
    tlc_wpp_dict.update(n_dict)

    fname = os.path.join('..', 'refdat_namesets', 'tlc_wpp_groups.json')
    with open(fname) as fid01:
        n_dict = json.load(fid01)
    tlc_wpp_dict.update(n_dict)

    # Parse CSVs
    wppf1 = 'WPP2022_POP_F02_1_POPULATION_BY_AGE_BOTH_SEXES_ESTIMATES.csv'
    with open(wppf1, errors='ignore') as fid01:
        flines_rev = [val.strip().split(',') for val in fid01.readlines()]

    wppf2 = 'WPP2022_POP_F02_1_POPULATION_BY_AGE_BOTH_SEXES_MEDIUM_VARIANT.csv'
    with open(wppf2, errors='ignore') as fid01:
        flines_fwd = [val.strip().split(',') for val in fid01.readlines()]

    # Construct output data structure
    rng = range(11, 32)
    pop_dat = np.zeros((0, 22), dtype=int)

    # Add values from retrospective estimates
    for rval in flines_rev:
        if (reprule(rval[2]) == tlc_wpp_dict[TLC]):
            year_val = int(rval[10])
            if (year_val % 5):
                continue
            bpop = [int(1000*float(rval[idx].replace(' ', ''))) for idx in rng]
            pop_dat = np.vstack((pop_dat, np.array([year_val]+bpop)))

    # Add values from forward projections
    for rval in flines_fwd:
        if (reprule(rval[2]) == tlc_wpp_dict[TLC]):
            year_val = int(rval[10])
            if (year_val % 5):
                continue
            if (year_val == pop_dat[-1, 0]):
                continue
            bpop = [int(1000*float(rval[idx].replace(' ', ''))) for idx in rng]
            pop_dat = np.vstack((pop_dat, np.array([year_val]+bpop)))

    # Write data files
    ofile_name = 'pop_dat_{:s}.csv'.format(TLC)
    np.savetxt(ofile_name, pop_dat.T, fmt='%d', delimiter=',')

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_pop_dat()

# *****************************************************************************
