# *****************************************************************************
#
# *****************************************************************************

import global_data as gdata

# *****************************************************************************


def application(timestep):

    TIME_VAL = float(timestep)

    # Example interface for in-processing;
    if (gdata.first_call_bool):
        gdata.first_call_bool = False

        msg_str = 'Hello and goodbye from in-process at time {:.1f}'
        print(msg_str.format(TIME_VAL))

    return None

# *****************************************************************************
