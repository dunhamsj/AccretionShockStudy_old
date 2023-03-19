#!/usr/bin/env python3

import numpy as np
from os.path import isdir

from ComputePowerInLegendreModes import ComputePowerInLegendreModes
from FitPowerToModel import FitPowerToModel
from computeTimeScales import ComputeTimeScales

R    = np.array( [ 'NR', 'GR' ], str )
M    = np.array( [ '2.8' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '6.00e1', '7.50e1', '9.00e1' ], str )
R    = np.array( [ 'NR' ], str )
Rs   = np.array( [ '7.50e1'], str )

R    = np.array( [ 'GR' ], str )
M    = np.array( [ '1.4' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '180' ], str )

for r in range( R.shape[0] ):
    for m in range( M.shape[0] ):
        for mdot in range( Mdot.shape[0] ):
            for rs in range( Rs.shape[0] ):

                #ID = '{:}2D_M{:}_Mdot{:}_Rs{:}'.format \
                #     ( R[r], M[m], Mdot[mdot], Rs[rs] )
                ID = '{:}2D_M{:}_Rpns040_Rs{:}_Mdot{:}'.format \
                     ( R[r], M[m], Rs[rs], Mdot[mdot] )

                #plotFileDirectory \
                #  = '/lump/data/accretionShockStudy/{:}/'.format( ID )
                plotFileDirectory \
                  = '/lump/data/accretionShockStudy/newRuns/newProductionRuns/{:}/'.format( ID )

                if not isdir( plotFileDirectory ):
                    print( '\n{:} does not exist. Skipping.\n' \
                           .format( plotFileDirectory ) )
                    continue

                #plotFileBaseName = '{:}.plt_'.format( ID )
                plotFileBaseName = '{:}.plt'.format( ID )

                dataFileName = '.{:}_LegendrePowerSpectrum.dat'.format( ID )

                ComputePowerInLegendreModes \
                  ( plotFileDirectory, plotFileBaseName, dataFileName, \
                    'DivV2', 0.8, 0.9, np.float64( Rs[rs] ), verbose = True )

                t, P0, P1, P2, P3, P4 = np.loadtxt( dataFileName )

                LogF  = np.log( 1.0e14 )
                tauR  = 200.0
                delta = 0.0

                rInner = 4.00e1
                rOuter = np.float64( Rs[rs] )
                tAd, tAc \
                  = ComputeTimeScales \
                      ( plotFileDirectory+plotFileBaseName+'00000000', \
                        rInner, rOuter, R[r] )

                T_SASI = tAd + tAc
                tF0 = 1.0
                tF1 = 120.0

                if M[m] == '1.4':
                    if Rs[rs] == '120':
                        T_SASI = 25.0
                        tF0    = 63.5
                        tF1    = 140.0
                    elif Rs[rs] == '150':
                        T_SASI = 35.0
                        tF0    = 25.0
                        tF1    = 140.0
                    elif Rs[rs] == '180':
                        T_SASI = 55.0
                        tF0    = 35.0
                        tF1    = 140.0
                elif M[m] == '2.0':
                    if Rs[rs] == '120':
                        T_SASI = 20.0
                        tF0    = 1.0
                        tF1    = 150.0
                    elif Rs[rs] == '150':
                        T_SASI = 30.0
                        tF0    = 20.0
                        tF1    = 140.0
                    elif Rs[rs] == '180':
                        T_SASI = 50.0
                        tF0    = 1.0
                        tF1    = 150.0
                elif M[m] == '2.8':
                    if Rs[rs] == '120':
                        T_SASI = 20.0
                        tF0    = 15.0
                        tF1    = 150.0
                    elif Rs[rs] == '150':
                        T_SASI = 30.0
                        tF0    = 55.0
                        tF1    = 150.0
                    elif Rs[rs] == '180':
                        T_SASI = 40.0
                        tF0    = 5.0
                        tF1    = 150.0
                    if R[r] == 'NR':
                        if Rs[rs] == '6.00e1':
                            tF0    = 1.0
                            tF1    = 50.0
                        elif Rs[rs] == '7.50e1':
                            tF0    = 1.0
                            tF1    = 90.0

                omegaR = 2.0 * np.pi / tauR
                omegaI = 2.0 * np.pi / T_SASI

                InitialGuess \
                  = np.array( [ LogF, omegaR, omegaI, delta ], np.float64 )

                dataFileName = '.{:}_Fit.dat'.format( ID )

                FitPowerToModel( tF0, tF1, t, P1, InitialGuess, dataFileName )

import os
os.system( 'rm -rf __pycache__ ' )
