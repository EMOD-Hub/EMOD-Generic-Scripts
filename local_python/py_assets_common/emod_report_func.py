# *****************************************************************************
#
# *****************************************************************************


from emod_constants import RST_FILE

# *****************************************************************************


def report_strain(report_dict, every_timestep=False):

    rst_str = RST_FILE.split('.')[0]

    report_dict['Custom_Reports'][rst_str] = {'Enabled': 1,
                                              'Reports': list()}

    rep_dict = {'Report_Name': RST_FILE}
    if (every_timestep):
        rep_dict['Ouput_Every_Timestep'] = 1

    report_dict['Custom_Reports'][rst_str]['Reports'].append(rep_dict)

    return None

# *****************************************************************************
