# *****************************************************************************
#
# *****************************************************************************

import global_data as gdata

# *****************************************************************************


def application(timestep):

    # Example interface for in-processing;
    if (gdata.first_call_bool):
        gdata.first_call_bool = False

        timeval = float(timestep)
        msg_str = 'Hello and goodbye from in-process at time {:.1f}'
        print(msg_str.format(timeval))

    return None

# *****************************************************************************
