#!/usr/bin/env python3

import numpy as np
from os.path import isdir

from ComputePowerInLegendreModes import ComputePowerInLegendreModes
from FitDataToModel import FitDataToModel
from computeTimeScales import ComputeTimeScales

def leg():

    # Early-stage
    R    = np.array( [ 'NR', 'GR' ], str )
    M    = np.array( [ '1.4', '2.0', '2.8' ], str )
    Mdot = np.array( [ '0.3' ], str )
    Rs   = np.array( [ '120', '150', '180' ], str )
    RPNS = 4.0e1
    suffix = ''

    # Late-stage
    R    = np.array( [ 'NR', 'GR' ], str )
    M    = np.array( [ '2.8' ], str )
    Mdot = np.array( [ '0.3' ], str )
    Rs   = np.array( [ '6.00e1', '7.50e1', '9.00e1' ], str )
    RPNS = 2.0e1
    suffix = '_RPNS2.00e1'

    R    = np.array( [ 'NR' ], str )
    M    = np.array( [ '2.8' ], str )
    Mdot = np.array( [ '0.3' ], str )
    Rs   = np.array( [ '6.00e1' ], str )

    for r in range( R.shape[0] ):
        for m in range( M.shape[0] ):
            for mdot in range( Mdot.shape[0] ):
                for rs in range( Rs.shape[0] ):

                    ID = '{:}2D_M{:}_Mdot{:}_Rs{:}{:}'.format \
                         ( R[r], M[m], Mdot[mdot], Rs[rs], suffix )

                    plotFileDirectory \
                      = '/lump/data/accretionShockStudy/{:}/'.format( ID )

                    if not isdir( plotFileDirectory ): continue

                    plotFileBaseName = '{:}.plt'.format( ID )

                    dataFileName = '.{:}_LegendrePowerSpectrum.dat'.format( ID )

                    if R[r] == 'NR':
                        f = 'LateralMomentumFluxInRadialDirectionNR'
                    else:
                        f = 'LateralMomentumFluxInRadialDirectionGR'

                    ComputePowerInLegendreModes \
                      ( plotFileDirectory, plotFileBaseName, dataFileName, \
                        f, 0.8, 0.9, np.float64( Rs[rs] ), verbose = True )

                    t, P0, P1, P2, P3, P4 = np.loadtxt( dataFileName )

                    LogF0 = np.log( P1[0] )
                    tau   = 5.0 # ms
                    delta = 0.0

                    rInner = RPNS
                    rOuter = np.float64( Rs[rs] )
                    tAd, tAc \
                      = ComputeTimeScales \
                          ( plotFileDirectory+plotFileBaseName+'00000000', \
                            rInner, rOuter, R[r] )
                    T   = tAd + tAc

                    tF0 = t[3]
                    tF1 = t[np.where( t < 5.0e1 )[0][-1]]

                    omegaR = 1.0 / tau
                    omegaI = 2.0 * np.pi / T

                    InitialGuess \
                      = np.array( [ LogF0, omegaR, omegaI, delta ], np.float64 )

                    dataFileName = '.{:}_Fit.dat'.format( ID )

                    dfn = '.{:}_LegendrePowerSpectrum.dat2.dat'.format( ID )
                    t, P1 = np.loadtxt( dfn )
                    InitialGuess \
                      = np.array( [ abs(P1[1]), omegaR, omegaI, delta ], \
                                  np.float64 )
                    FitDataToModel \
                      ( tF0, tF1, t, P1, InitialGuess, dataFileName )

leg()

import os
os.system( 'rm -rf __pycache__ ' )
