#!/usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit

def FitPowerToModel( t0, t1, t, P1, InitialGuess, dataFileName ):

    ind = np.where( ( t >= t0 ) & ( t <= t1 ) )[0]

    tFit = t[ind] - t0

    beta, pcov = curve_fit( FittingFunction, tFit, np.log( P1[ind] ), \
                            p0 = InitialGuess, jac = Jacobian )

    perr = np.sqrt( np.diag( pcov ) )

    np.savetxt \
      ( dataFileName, \
        [ t0, t1, \
          beta[0], beta[1], beta[2], beta[3], \
          perr[0], perr[1], perr[2], perr[3] ], \
          header = \
't0, t1, LogF1, omegaR, omegaI, delta, dLogF1, domegaR, domegaI, ddelta' )

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

