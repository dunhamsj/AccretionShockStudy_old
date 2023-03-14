#!/usr/bin/env python3

import numpy as np
from os.path import isdir
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from FitPowerToModel import FittingFunction
from computeTimeScales import ComputeTimeScales

R    = np.array( [ 'NR', 'NR', 'NR', 'NR' ], str )
M    = np.array( [ '1.4' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '180' ], str )

arrShape = (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0])

ID = np.empty( arrShape, object )

t  = np.empty( arrShape, object )
P0 = np.empty( arrShape, object )
P1 = np.empty( arrShape, object )
P2 = np.empty( arrShape, object )
P3 = np.empty( arrShape, object )
P4 = np.empty( arrShape, object )

t0         = np.empty( arrShape, object )
t1         = np.empty( arrShape, object )
LogF       = np.empty( arrShape, object )
omegaR     = np.empty( arrShape, object )
omegaI     = np.empty( arrShape, object )
delta      = np.empty( arrShape, object )
omegaR_err = np.empty( arrShape, object )
omegaI_err = np.empty( arrShape, object )

T_SASI     = np.empty( arrShape, object )

for r in range( R.shape[0] ):
    for m in range( M.shape[0] ):
        for mdot in range( Mdot.shape[0] ):
            for rs in range( Rs.shape[0] ):

                if r == 0:

                    ID[r,m,mdot,rs] \
                      = '{:}2D_M{:}_Mdot{:}_Rs{:}'.format \
                         ( R[r], M[m], Mdot[mdot], Rs[rs] )

                    plotFileDirectory \
                      = '/lump/data/accretionShockStudy/{:}/'.format \
                        ( ID[r,m,mdot,rs] )

                else:

                    ID[r,m,mdot,rs] \
                      = '{:}2D_M{:}_Rpns040_Rs{:}_Mdot{:}'.format \
                         ( R[r], M[m], Rs[rs], Mdot[mdot] )

                    plotFileDirectory \
                      = '/lump/data/accretionShockStudy/newRuns/newProductionRuns/{:}/'.format \
                        ( ID[r,m,mdot,rs] )

                if not isdir( plotFileDirectory ):
                    print( '{:} does not exist. Skipping.' \
                           .format( plotFileDirectory ) )
                    continue

                if r == 0:
                    plotFileBaseName = '{:}.plt_'.format( ID[r,m,mdot,rs] )
                else:
                    plotFileBaseName = '{:}.plt'.format( ID[r,m,mdot,rs] )

                rInner = 4.00e1
                rOuter = np.float64( Rs[rs] )
                tAd, tAc \
                  = ComputeTimeScales \
                      ( plotFileDirectory+plotFileBaseName+'00000000', \
                        rInner, rOuter, R[r] )

                T_SASI[r,m,mdot,rs] = tAd + tAc

                if r == 0:
                    dataFileName \
                      = '.{:}_LegendrePowerSpectrum_original.dat'.format( ID[r,m,mdot,rs] )
                if r == 1:
                    dataFileName \
                      = '.{:}_LegendrePowerSpectrum_PA1.0e-05.dat'.format( ID[r,m,mdot,rs] )
                if r == 2:
                    dataFileName \
                      = '.{:}_LegendrePowerSpectrum_oldPert.dat'.format( ID[r,m,mdot,rs] )
                if r == 3:
                    dataFileName \
                      = '.{:}_LegendrePowerSpectrum_origPert.dat'.format( ID[r,m,mdot,rs] )

                t [r,m,mdot,rs], \
                P0[r,m,mdot,rs], \
                P1[r,m,mdot,rs], \
                P2[r,m,mdot,rs], \
                P3[r,m,mdot,rs], \
                P4[r,m,mdot,rs] \
                  = np.loadtxt( dataFileName )

                dataFileName = '.{:}_Fit.dat'.format( ID[r,m,mdot,rs] )

                t0        [r,m,mdot,rs], \
                t1        [r,m,mdot,rs], \
                LogF      [r,m,mdot,rs], \
                omegaR    [r,m,mdot,rs], \
                omegaI    [r,m,mdot,rs], \
                delta     [r,m,mdot,rs], \
                dummy0, \
                omegaR_err[r,m,mdot,rs], \
                omegaI_err[r,m,mdot,rs], \
                dummy1 \
                  = np.loadtxt( dataFileName )

fig, ax = plt.subplots( 1, 1 )#, figsize = (12,9) )

mdot = 0

lab = [ 'Original', \
        r'$\Delta p/p=10^{{-5}}$', \
        'Old Perturbation Method (xH = 270 km)', \
        'Old Perturbation Method (xH = 360 km)' ]

for r in range( R.shape[0] ):
    for m in range( M.shape[0] ):
        for rs in range( Rs.shape[0] ):

            tt  = t [r,m,mdot,rs]
            P11 = P1[r,m,mdot,rs]

            ax.plot( tt, P11, '-', label = lab[r] )

ax.grid()
ax.set_yscale( 'log' )
ax.set_xlim( 0, 1.0e2 )
ax.set_title \
  ( r'$\texttt{{NR2D_M{:}_Rpns040_Rs180_Mdot0.3}}$' \
    .format( M[0] ), fontsize = 15 )

ax.legend()
fig.supxlabel( 'Time [ms]' )
fig.supylabel( r'$H_{1}$ [cgs]' )

plt.savefig( '/home/kkadoogan/fig.PowerInLegendreMode.png', \
             dpi = 300 )
#plt.show()

import os
os.system( 'rm -rf __pycache__ ' )
