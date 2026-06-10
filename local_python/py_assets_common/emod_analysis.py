# *****************************************************************************
#
# *****************************************************************************

import numpy as np

# *****************************************************************************

EPS = np.finfo(float).eps

# *****************************************************************************


def norpois_vec(yobs: list[float], ysim: list[float],
                mval: float = 0.0, bval: float = 1.0):
    '''
    Poisson-based estimation of log-liklihood for a timeseries. The vector of
    observed values is used as the occurances (k) and the vector of simulated
    values is used as the rates (lam). Stirling's approximation is applied to
    the probability:
        ln P = k*ln(lam/k) - lam - k - 0.5*ln(2*pi*k)

    A small positive value is added to the simulated rates so they are never
    exactly zero. The adjusted rates are then modified by a scale factor
    that represents the case to infection ratio. Scale factor adjustements
    have an `m*x + b` form that allows for linear changes over time. The values
    in the timeseries are assumed to be equally spaced, so the `x` value in
    the scale factor is the index of the value in the vector.
        lam = (ysim + delt)*(mval*x + bval)
        k = yobs

    First and second derivatives of the `mval` and `bval` coefficients in the
    scale factor are also calculated. Derivatives of the `bval` coefficient are
    with respect to `ln(bval)` to ensure positivity.

    Args:

    Returns:

    '''

    Ltot = 0.0
    Gtot = np.zeros((2,1))
    Htot = np.zeros((2,2))
    delt = 0.1

    for k1 in range(len(yobs)):
        yobs0 = float(yobs[k1])
        ysim0 = float(ysim[k1]) + delt
        rfac = mval*k1 + bval

        if (rfac <= 0.0):
            print(k1, mval, bval)
            raise ValueError("Observation rate is negative.")

        if (yobs0 > 0):
            Ltot += yobs0*np.log((ysim0*rfac)/yobs0) - \
                    ysim0*rfac + yobs0 - 0.5*np.log(2.0*np.pi*yobs0)

            Gtot[0, 0] += yobs0*k1/rfac - ysim0*k1
            Gtot[1, 0] += yobs0*bval/rfac - ysim0*bval

            Htot[0, 0] += -yobs0*k1*k1/rfac/rfac
            Htot[1, 0] += -yobs0*bval*k1/rfac/rfac
            Htot[0, 1] += -yobs0*bval*k1/rfac/rfac
            Htot[1, 1] += yobs0*bval/rfac - yobs0*bval*bval/rfac/rfac - \
                          ysim0*bval

        elif (yobs0 == 0):
            Ltot += -ysim0*rfac

            Gtot[0, 0] += -ysim0*k1
            Gtot[1, 0] += -ysim0*bval

            Htot[1, 1] += -ysim0*bval

    dval = np.linalg.solve(Htot, Gtot)

    return (Ltot, dval)

# *****************************************************************************


def norpois_opt(yobs: list[float], ysim: list[float]):
    '''
    '''

    nobs = len(yobs)
    nsim = len(ysim)

    if (nsim < nobs):
        raise ValueError("Simulation vector is too short.")

    sval = np.zeros((2,1))

    while (True):
        (Ltot, dval) = norpois_vec(yobs, ysim, sval[0, 0],
                                    np.exp(sval[1, 0]))

        # Convergence check
        dmag = np.max(np.abs(dval))
        if (dmag < 1.0e-4):
            break

        sval_new = sval - dval
        rfac_new = sval_new[0, 0]*nobs + np.exp(sval_new[1, 0])
        while (rfac_new <= 0):
            dval = dval / 3.0
            sval_new = sval - dval
            rfac_new = sval_new[0, 0]*nobs + np.exp(sval_new[1, 0])

        sval = sval_new
        print(dval)
        print(sval)
        print()

    # Adjust intercept to linear scale
    sval[1, 0] = np.exp(sval[1, 0])

    return (Ltot, (sval[0,0], sval[1, 0]))

# *****************************************************************************


def gauss_vec(dmu2, sig2, yscal=1.0):

    llik = -0.5*np.sum(dmu2/(sig2*yscal*yscal) +
                       np.log(2*np.pi*sig2*yscal*yscal))
    G = 1.0*np.sum(dmu2/(sig2*yscal*yscal)-1.0)
    H = -2.0*np.sum(dmu2/(sig2*yscal*yscal))

    return (llik, G, H)

# *****************************************************************************


def gauss_opt(yobs, ysim):

    k = np.array(yobs[0])
    N = np.array(yobs[1])
    mu = np.array(ysim[0])/(np.array(ysim[1])+EPS)
    x = k/(N+EPS)
    sig2 = 1.0/N[N > 0]
    dmu2 = np.power(x-mu, 2.0)[N > 0]

    lyscal = 0.0

    while (True):
        (llik, G, H) = gauss_vec(dmu2, sig2, np.exp(lyscal))
        sstptot = G/(H+EPS)

        # Step size control
        if (abs(sstptot) > 5.0):
            sstptot = np.sign(sstptot)*5.0

        lyscal = lyscal - sstptot

        if (abs(sstptot) < 1.0e-4):
            break

    return (llik, np.exp(lyscal))

# *****************************************************************************


def binom_vec(yobs, ysim):

    lliktot = 0

    pVec = ysim
    kVec = yobs[0]
    NVec = yobs[1]

    for k1 in range(len(kVec)):
        p = min(1-pVec[k1], 1-1e-7)
        k = kVec[k1]
        N = NVec[k1]

        if (N == 0):
            continue

        if ((p == 0 and k > 0) or (p == 1 and k < N)):
            raise Exception('Bad binomial data')

        if (k == 0):
            lliktot = lliktot + N*np.log(1-p)
            continue

        if (k == N):
            lliktot = lliktot + k*np.log(p)
            continue

        lliktot += k*np.log(p+1e-10) + (N-k)*np.log(1-p+1e-10) + \
            N*np.log(N/(N-k)) - k*np.log(k/(N-k)) + \
            0.5*np.log(N/k/2/np.pi/(N-k))

    return lliktot

# *****************************************************************************


def multinom_vec(yobs, ysim):

    lliktot = 0.0
    mlam = 0.1
    yobstot = 0.0
    ysimtot = 0.0

    for k1 in range(len(yobs)):
        yobsval = float(yobs[k1])
        ysimval = float(ysim[k1])

        if (yobsval >= 0):
            lliktot = lliktot + yobsval*np.log((ysimval + mlam)/yobsval) \
                              - 0.5*np.log(2.0*np.pi*yobsval)
            yobstot = yobstot + yobsval
            ysimtot = ysimtot + ysimval + mlam

    if (yobstot > 0):
        lliktot = lliktot + yobstot*np.log(yobstot/ysimtot) + \
                              0.5*np.log(2.0*np.pi*yobstot)

    return lliktot

# *****************************************************************************
