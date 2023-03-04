#!/usr/bin/env python3

import numpy as np
from os.path import isdir

from ComputePowerInLegendreModes import ComputePowerInLegendreModes
from FitPowerToModel import FitPowerToModel
from computeTimeScales import ComputeTimeScales

R    = 'NR'
M    = '2.8'
Mdot = '0.3'
Rs   = '120'
nX2  = np.array( [ '064' ], str )
nX2  = np.array( [ '128', '256' ], str )

for N in range( nX2.shape[0] ):

    ID = '{:}2D_M{:}_Mdot{:}_Rs{:}_nX768x{:}'.format \
         ( R, M, Mdot, Rs, nX2[N] )

    plotFileDirectory \
      = '/lump/data/accretionShockStudy/angularResolution/{:}/'.format( ID )

    if not isdir( plotFileDirectory ): continue

    plotFileBaseName = '{:}.plt_'.format( ID )

    dataFileName = '.{:}_LegendrePowerSpectrum.dat'.format( ID )

    ComputePowerInLegendreModes \
      ( plotFileDirectory, plotFileBaseName, dataFileName, \
        'DivV2', 0.8, 0.9, np.float64( Rs ), verbose = True )

    t, P0, P1, P2, P3, P4 = np.loadtxt( dataFileName )

    LogF  = np.log( 1.0e14 )
    tauR  = 200.0
    delta = 0.0

    rInner = 4.00e1
    rOuter = np.float64( Rs )
    tAd, tAc \
      = ComputeTimeScales \
          ( plotFileDirectory+plotFileBaseName+'00000000', \
            rInner, rOuter, R )

    T_SASI = tAd + tAc
    tF0 = 1.0
    tF1 = 100.0

    omegaR = 2.0 * np.pi / tauR
    omegaI = 2.0 * np.pi / T_SASI

    InitialGuess \
      = np.array( [ LogF, omegaR, omegaI, delta ], np.float64 )

    dataFileName = '.{:}_Fit.dat'.format( ID )

    FitPowerToModel( tF0, tF1, t, P1, InitialGuess, dataFileName )

import os
os.system( 'rm -rf __pycache__ ' )
