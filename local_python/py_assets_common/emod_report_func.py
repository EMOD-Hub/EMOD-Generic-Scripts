# *****************************************************************************
#
# *****************************************************************************

import json

from emod_api import schema_to_class as s2c

from emod_constants import REPORTS_FILE, REPORTS_KEY, SPATH, RST_FILE

# *****************************************************************************


def build_file():

    json_set = {REPORTS_KEY: list()}

    return json_set

# *****************************************************************************


def write_file(json_set):

    for rep_obj in json_set[REPORTS_KEY]:
        rep_obj.finalize()

    with open(REPORTS_FILE, 'w') as fid01:
        json.dump(json_set, fid01, sort_keys=True, indent=4)

    return None

# *****************************************************************************


def report_strain(json_set, time_start=0.0, every_timestep=False):

    rst_class = RST_FILE.split('.')[0]
    rep_dict = s2c.get_class_with_defaults(rst_class, SPATH)

    rep_dict.Time_Start = time_start

    if (every_timestep):
        rep_dict.Ouput_Every_Timestep = 1

    json_set[REPORTS_KEY].append(rep_dict)

    return None

# *****************************************************************************
