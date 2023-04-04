#!/usr/bin/env python3

import numpy as np
from os.path import isfile, isdir
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from FitPowerToModel import FittingFunction
from computeTimeScales import ComputeTimeScales

R    = np.array( [ 'NR' ], str )
M    = np.array( [ '1.4' ], str )
Rs   = np.array( [ '1.20e2' ], str )

arrShape = (R.shape[0],M.shape[0],Rs.shape[0])

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
        for rs in range( Rs.shape[0] ):

            ID[r,m,rs] \
              = '{:}2D_M{:}_Rpns040_Rs{:}'.format \
                 ( R[r], M[m], Rs[rs] )

            plotFileDirectory \
              = '/lump/data/accretionShockStudy/newData/2D/{:}/'.format \
                ( ID[r,m,rs] )

            if not isdir( plotFileDirectory ):
                print( '{:} does not exist. Skipping.' \
                       .format( plotFileDirectory ) )
                continue

            plotFileBaseName = '{:}.plt'.format( ID[r,m,rs] )
            rInner = 4.00e1
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

fig, ax = plt.subplots( 1, 1 )#, figsize = (12,9) )

mdot = 0
for r in range( R.shape[0] ):
    for m in range( M.shape[0] ):
        for rs in range( Rs.shape[0] ):

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

            ind = np.where( ( tt >= t0t ) & ( tt <= t1t ) )[0]

            tF = tt[ind]

            logFt   = LogF  [r,m,rs]
            omegaRt = omegaR[r,m,rs]
            omegaIt = omegaI[r,m,rs]
            deltat  = delta [r,m,rs]

#            F = FittingFunction \
#                 ( tF - tF[0], logFt, \
#                   omegaRt, omegaIt, deltat )

            tau = T_SASI[0,m,rs]

            ind = np.where( ( tt < 601.0 ) & ( tt >= 0.0 ) )[0]

            if r == 1:
                ax.plot \
                  ( tt[ind]/tau, P1t[ind], '--', color = color[rs], \
                    label = 'Rs={:} km (GR)'.format( Rs[rs] ) )

            else:
                ax.plot \
                  ( tt[ind]/tau, P1t[ind], '-', color = color[rs], \
                    label = 'Rs={:} km'.format( Rs[rs] ) )
ax.grid()
ax.set_yscale( 'log' )
#ax.set_xlim( 0, 12 )
#xticks = np.linspace( 0, 12, 13, dtype = np.int64 )
#xticklabels = [ str( i ) for i in xticks ]
#ax.set_xticks( xticks )
#ax.set_xticklabels( xticklabels )
##ax.axvline( 150/T_SASI[0,0,0,0], label = r'$t=150\,\mathrm{ms}$' )

ax.set_title( r'$\texttt{{2D_M{:}_Rpns040}}$'.format( M[0] ), \
              fontsize = 15 )

ax.legend()
fig.supxlabel( r'$t/T_{\mathrm{SASI,NR}}$' )
#fig.supxlabel( 'Time [ms]' )
fig.supylabel( r'$H_{1}$ [cgs]' )
plt.subplots_adjust( hspace = 0.0, wspace = 0.3 )

#plt.savefig( '/home/kkadoogan/fig.PowerInLegendreMode.png', dpi = 300 )
plt.show()

import os
os.system( 'rm -rf __pycache__ ' )
