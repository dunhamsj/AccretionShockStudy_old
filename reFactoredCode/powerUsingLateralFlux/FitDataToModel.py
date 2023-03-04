#!/usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit

def FitDataToModel( t0, t1, t, Data, InitialGuess, dataFileName ):

    ind = np.where( ( t >= t0 ) & ( t <= t1 ) )[0]

    tFit = t[ind] - t0

    beta, pcov = curve_fit( FittingFunction, tFit, Data[ind], \
                            p0 = InitialGuess, jac = Jacobian )

    perr = np.sqrt( np.diag( pcov ) )

    np.savetxt \
      ( dataFileName, \
        [ t0, t1, \
          beta[0], beta[1], beta[2], beta[3], \
          perr[0], perr[1], perr[2], perr[3] ], \
          header = \
't0, t1, logF1, omegaR, omegaI, delta, dlogF1, domegaR, domegaI, ddelta' )

def FittingFunction( t, logF1, omegaR, omegaI, delta ):

    #  return logF1 + omegaR * t \
    #           + np.log( np.sin( omegaI * t + delta ) )
  return np.exp( logF1 ) * np.exp( omegaR * t ) * np.sin( omegaI * t + delta )

def Jacobian( t, logF1, omegaR, omegaI, delta ):

  # Jacobian of fitting function

  J = np.empty( (t.shape[0],4), np.float64 )

  Phase = omegaI * t + delta

  A = np.exp( logF1 )
  B = np.exp( omegaR * t )
  C = np.sin( Phase )
  D = np.cos( Phase )

  J[:,0] = B * C
  J[:,1] = A * B * C * t
  J[:,2] = A * B * D * t
  J[:,3] = A * B * D
#  J[:,0] = 1.0
#
#  J[:,1] = t
#
#  J[:,2] = 1.0 / np.sin( Phase ) * np.cos( Phase ) * t
#
#  J[:,3] = 1.0 / np.sin( Phase ) * np.cos( Phase )

  return J

