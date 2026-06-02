# *****************************************************************************
#
# *****************************************************************************

import numpy as np

# *****************************************************************************

    """
    Helper function to calculate approximate gravity model coefficients for a
    list of gravity model exponents. Approximate fit using rough empirical
    correlation. Ought to be tuned with data; should be adequate as an initial
    guess.

    Args:
        exp_vals (list[float]): List of gravity model exponent values.

        ln_offset (list[float]): List of offsets for each exponent. Adjustment
            assumes log-scale offset.

    Returns:
        (list[float]): Estimated gravity model coffefficient for each gravity
            model exponent provided as an input argument.
     """


def grav_coeff_guess(exp_vals: list[float], ln_offset: list[float] = None):

    if (np.min(exp_vals) < 0.0 or np.max(exp_vals) > 8.0):
        raise Exception('Network exponent out of range.')

    x_ref = np.array([0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 7, 8])
    y_ref = np.array([-2.794, -1.298, 0.155, 1.528, 2.797, 6.924,
                      9.774, 12.22, 14.44, 16.56, 18.65, 20.75])

    ln_coeffs = np.interp(exp_vals, x_ref, y_ref) + np.array(ln_offset)


    ni_coeff = np.exp(ln_coeffs).to_list()

    return (np.exp(ln_coeffs)).to_list()

# *****************************************************************************
