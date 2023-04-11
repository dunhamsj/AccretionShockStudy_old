#!/usr/bin/env python3

import numpy as np
from os.path import isfile, isdir
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from FitPowerToModel import FittingFunction
from computeTimeScales import ComputeTimeScales

R    = np.array( [ 'NR', 'GR' ], str )
M    = np.array( [ '1.4' ], str )
Rs   = np.array( [ '1.20e2', '1.50e2', '1.75e2' ], str )
Rpns = np.array( [ '040' ], str )

arrShape = (R.shape[0],M.shape[0],Rs.shape[0])

ID = np.empty( arrShape, object )

indd = np.empty( arrShape, np.int64 )
indd[0,0,0] = 777
indd[0,0,1] = 1330
indd[0,0,2] = 1839
indd[1] = np.copy( indd[0] )

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
        for rs in range( Rs.shape[0] ):

            ID[r,m,rs] \
              = '{:}2D_M{:}_Rpns{:}_Rs{:}'.format \
                 ( R[r], M[m], Rpns[0], Rs[rs] )

            plotFileDirectory \
              = '/lump/data/accretionShockStudy/newData/2D/{:}/'.format \
                ( ID[r,m,rs] )

            if not isdir( plotFileDirectory ):
                print( '{:} does not exist. Skipping.' \
                       .format( plotFileDirectory ) )
                continue

            plotFileBaseName = '{:}.plt'.format( ID[r,m,rs] )
            rInner = np.float64( Rpns[0] )
            rOuter = np.float64( Rs[rs] )
            tAd = 0
            tAc = 0
            tAd, tAc \
              = ComputeTimeScales \
                  ( plotFileDirectory+plotFileBaseName+'00000000', \
                    rInner, rOuter, R[r] )

            T_SASI[r,m,rs] = tAd + tAc

            dataFileName \
              = '.{:}_LegendrePowerSpectrum.dat'.format( ID[r,m,rs] )

            if not isfile( dataFileName ):
                print( '{:} does not exist. Skipping.' \
                       .format( dataFileName ) )
                continue

            t [r,m,rs], \
            P0[r,m,rs], \
            P1[r,m,rs], \
            P2[r,m,rs], \
            P3[r,m,rs], \
            P4[r,m,rs] \
              = np.loadtxt( dataFileName )

            # Read in fit data

            f = open( dataFileName )

            dum = f.readline()

            s = f.readline(); ind = s.find( '#' )+1
            tmp \
              = np.array( list( map( np.float64, s[ind:].split() ) ), \
                          np.float64 )
            t0    [r,m,rs] = tmp[0]
            t1    [r,m,rs] = tmp[1]
            LogF  [r,m,rs] = tmp[2]
            omegaR[r,m,rs] = tmp[3]
            omegaI[r,m,rs] = tmp[4]
            delta [r,m,rs] = tmp[5]

            dum = f.readline()

            s = f.readline(); ind = s.find( '#' )+1
            tmp \
              = np.array( list( map( np.float64, s[ind:].split() ) ), \
                          np.float64 )
            omegaR_err[r,m,rs] = tmp[1]
            omegaI_err[r,m,rs] = tmp[2]

            f.close()

# colorblind-friendly palette: https://gist.github.com/thriveth/8560036
color = ['#377eb8', '#ff7f00', '#4daf4a', \
         '#f781bf', '#a65628', '#984ea3', \
         '#999999', '#e41a1c', '#dede00']

fig, axs = plt.subplots( R.shape[0], Rs.shape[0], figsize = (16,9) )

mdot = 0
for r in range( R.shape[0] ):
    for rs in range( Rs.shape[0] ):
        for m in range( M.shape[0] ):

            dataFileName \
              = '.{:}_LegendrePowerSpectrum.dat'.format( ID[r,m,rs] )

            if not isfile( dataFileName ):
                print( '{:} does not exist. Skipping.' \
                       .format( dataFileName ) )
                continue

            tt  = t [r,m,rs]
            P1t = P1[r,m,rs]
            t0t = t0[r,m,rs]
            t1t = t1[r,m,rs]

            tt  = np.copy( tt [0:indd[r,m,rs]] )
            P1t = np.copy( P1t[0:indd[r,m,rs]] )

            ind = np.where( ( tt >= t0t ) & ( tt <= t1t ) )[0]

            tF = tt[ind]

            logFt   = LogF  [r,m,rs]
            omegaRt = omegaR[r,m,rs]
            omegaIt = omegaI[r,m,rs]
            deltat  = delta [r,m,rs]

            F = np.exp( FittingFunction \
                          ( tF - tF[0], logFt, \
                            omegaRt, omegaIt, deltat ) )

            tau = T_SASI[r,m,rs]

            axs[r,rs].plot( tt/tau, P1t, '-', color = color[rs] )
            axs[r,rs].plot( tF/tau, F  , '-', color = 'k' )

        axs[r,rs].grid()
        axs[r,rs].set_yscale( 'log' )
        axs[r,rs].set_xlim( -0.5, 10.5 )
        axs[r,rs].set_ylim( 1.0e10, 1.0e25 )
        axs[r,rs].set_title( r'$\texttt{{{:}}}$'.format( ID[r,m,rs] ), \
                             fontsize = 15 )

fig.supxlabel( r'$t/T_{\mathrm{SASI}}$', y = +0.05, fontsize = 15 )
#fig.supxlabel( 'Time [ms]' )
fig.supylabel( r'$H_{1}$ [cgs]', x = +0.075, fontsize = 15 )

plt.savefig( '/home/kkadoogan/fig.PowerInLegendreMode.png', dpi = 300 )
#plt.show()

import os
os.system( 'rm -rf __pycache__ ' )
