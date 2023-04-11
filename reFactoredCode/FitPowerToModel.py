#!/usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit

def FitPowerToModel( t0, t1, t, P1, InitialGuess ):

    ind = np.where( ( t >= t0 ) & ( t <= t1 ) )[0]

    tFit = t[ind] - t0

    LogF1_Lo  = 30.0
    LogF1_Hi  = 35.0
    omegaR_Lo = 1.0 / np.inf        # No growth of SASI
    omegaR_Hi = 1.0 / 1.0           # Power increase by factor of e in 1 ms
    omegaI_Lo = 2.0 * np.pi / 1.0e2 # Oscillatio period of 100 ms
    omegaI_Hi = 2.0 * np.pi / 1.0e0 # Oscillatio period of 1 ms
    delta_Lo  = 0.0
    delta_Hi  = np.pi

    bLo = [ LogF1_Lo, omegaR_Lo, omegaI_Lo, delta_Lo ]
    bHi = [ LogF1_Hi, omegaR_Hi, omegaI_Hi, delta_Hi ]
    beta, pcov \
      = curve_fit \
          ( FittingFunction, tFit, np.log( P1[ind] ), \
            p0 = InitialGuess, jac = Jacobian)#, bounds = ( bLo, bHi ) )

    perr = np.sqrt( np.diag( pcov ) )

    return beta, perr

def FittingFunction( t, logF1, omega_r, omega_i, delta ):

  # Modified fitting function
  # (log of Eq. (9) in Blondin & Mezzacappa, (2006))

  return logF1 + 2.0 * omega_r * t \
           + np.log( np.sin( omega_i * t + delta )**2 )

def Jacobian( t, logF1, omega_r, omega_i, delta ):

  # Jacobian of modified fitting function

  J = np.empty( (t.shape[0],4), np.float64 )

  ImPhase = omega_i * t + delta

  J[:,0] = 1.0

  J[:,1] = 2.0 * t

  J[:,2] = 2.0 * np.cos( ImPhase ) / np.sin( ImPhase ) * t

  J[:,3] = 2.0 * np.cos( ImPhase ) / np.sin( ImPhase )

  return J

