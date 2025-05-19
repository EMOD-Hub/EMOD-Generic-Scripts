# *****************************************************************************
#
# *****************************************************************************

import numpy as np

# *****************************************************************************

EPS = np.finfo(float).eps

# *****************************************************************************


def norpois_vec(yobs, ysim, yscal=1.0):

    lliktot = 0.0
    Gtot = 0.0
    Htot = 0.0
    mlam = 0.1

    for k1 in range(len(yobs)):
        yobsval = float(yobs[k1])
        ysimval = float(ysim[k1])
        llik = 0.0
        G = 0.0
        H = 0.0

        if (yobsval > 0):
            llik = yobsval*np.log((yscal*ysimval + mlam)/yobsval) \
                   - yscal*ysimval - mlam + yobsval \
                   - 0.5*np.log(2.0*np.pi*yobsval)
            G = yobsval - ysimval*yscal
            H = - ysimval*yscal
        elif (yobsval == 0):
            llik = -ysimval*yscal - mlam
            G = -ysimval*yscal
            H = -ysimval*yscal

        lliktot = lliktot + llik
        Gtot = Gtot + G
        Htot = Htot + H

    if (Htot != 0.0):
        sstptot = Gtot/Htot
    else:
        sstptot = 0.0

    return (lliktot, sstptot)

# *****************************************************************************


def norpois_opt(yobs, ysim):

    lyscal = 0.0

    while (True):
        (lliktot, sstptot) = norpois_vec(yobs, ysim, np.exp(lyscal))

        # Step size control
        if (abs(sstptot) > 5.0):
            sstptot = np.sign(sstptot)*5.0

        lyscal = lyscal - sstptot
        if (abs(sstptot) < 1.0e-4):
            break

    return (lliktot, np.exp(lyscal))

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
