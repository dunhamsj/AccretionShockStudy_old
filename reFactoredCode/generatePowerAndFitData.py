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

R    = np.array( [ 'NR', 'GR' ], str )
M    = np.array( [ '2.8' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '7.50e1' ], str )
Rpns = np.array( [ '020' ], str )

for r in range( R.shape[0] ):
    for m in range( M.shape[0] ):
        for mdot in range( Mdot.shape[0] ):
            for rs in range( Rs.shape[0] ):

                ID = '{:}2D_M{:}_Rpns{:}_Rs{:}'.format \
                     ( R[r], M[m], Rpns[0], Rs[rs] )

                plotFileDirectory \
                  = '/lump/data/accretionShockStudy/newData/2D/{:}/' \
                    .format( ID )

                if not isdir( plotFileDirectory ):
                    print( '\n{:} does not exist. Skipping.\n' \
                           .format( plotFileDirectory ) )
                    continue

                #plotFileBaseName = '{:}.plt_'.format( ID )
                plotFileBaseName = '{:}.plt'.format( ID )

                dataFileName = '.{:}_LegendrePowerSpectrum.dat'.format( ID )

                ComputePowerInLegendreModes \
                  ( plotFileDirectory, plotFileBaseName, dataFileName, \
                    'DivV2', 0.8, 0.9, np.float64( Rs[rs] ), \
                    fc = False, ow = False, verbose = True )

                t, P0, P1, P2, P3, P4 = np.loadtxt( dataFileName )

                LogF  = 31.0
                tauR  = 20.0
                delta = np.pi/4.0

                rInner = np.float64( Rpns[0] )
                rOuter = np.float64( Rs[rs] )
                tAd, tAc \
                  = ComputeTimeScales \
                      ( plotFileDirectory+plotFileBaseName+'00000000', \
                        rInner, rOuter, R[r] )

                T_SASI = tAd + tAc
                tF0 = 1.0 * T_SASI
                tF1 = 8.0 * T_SASI

                omegaR = 1.0 / tauR
                omegaI = 2.0 * np.pi / T_SASI

                InitialGuess \
                  = np.array( [ LogF, omegaR, omegaI, delta ], np.float64 )

                beta, perr \
                  = FitPowerToModel \
                      ( tF0, tF1, t, P1, InitialGuess )

                b = ''
                e = ''
                for i in range( len( beta ) ):
                    b += str( beta[i] ) + ' '
                for i in range( len( perr ) ):
                    e += str( perr[i] ) + ' '
                header = 'tF0, tF1, LogF1, omegaR, omegaI, delta\n' \
                         + str( tF0 ) + ' ' + str( tF1 ) + ' ' + b \
                         + '\ndLogF1, domegaR, domegaI, ddelta\n' + e \
                         + '\nTime [ms], P0 [cgs], P1 [cgs], ' \
                         + 'P2 [cgs], P3 [cgs], P4 [cgs]'

                t, P0, P1, P2, P3, P4 = np.loadtxt( dataFileName )

                Data = np.vstack( (t,P0,P1,P2,P3,P4) )
                np.savetxt( dataFileName, Data, header = header )

import os
os.system( 'rm -rf __pycache__ ' )
