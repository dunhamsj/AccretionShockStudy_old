#!/usr/bin/env python3

import numpy as np
from os.path import isdir
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from FitPowerToModel import FittingFunction
from computeTimeScales import ComputeTimeScales

R    = np.array( [ 'NR' ], str )
M    = np.array( [ '1.4' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '170', '175', '179', '180', '185', '190' ], str )

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

                #if r == 0:

                #    ID[r,m,mdot,rs] \
                #      = '{:}2D_M{:}_Mdot{:}_Rs{:}'.format \
                #         ( R[r], M[m], Mdot[mdot], Rs[rs] )

                #    plotFileDirectory \
                #      = '/lump/data/accretionShockStudy/{:}/'.format \
                #        ( ID[r,m,mdot,rs] )

                #else:

                #    ID[r,m,mdot,rs] \
                #      = '{:}2D_M{:}_Rpns040_Rs{:}_Mdot{:}'.format \
                #         ( R[r], M[m], Rs[rs], Mdot[mdot] )

                #    plotFileDirectory \
                #      = '/lump/data/accretionShockStudy/newRuns/{:}/'.format \
                #        ( ID[r,m,mdot,rs] )

                ID[r,m,mdot,rs] \
                  = '{:}2D_M{:}_Rpns040_Rs{:}_Mdot{:}'.format \
                     ( R[r], M[m], Rs[rs], Mdot[mdot] )

                plotFileDirectory \
                  = '/lump/data/accretionShockStudy/newRuns/newProductionRuns/{:}/'.format \
                    ( ID[r,m,mdot,rs] )

                #if not isdir( plotFileDirectory ):
                #    print( '{:} does not exist. Skipping.' \
                #           .format( plotFileDirectory ) )
                #    continue

                plotFileBaseName = '{:}.plt'.format( ID[r,m,mdot,rs] )
                rInner = 4.00e1
                rOuter = np.float64( Rs[rs] )
                tAd = 0
                tAc = 0
                #tAd, tAc \
                #  = ComputeTimeScales \
                #      ( plotFileDirectory+plotFileBaseName+'00000000', \
                #        rInner, rOuter, R[r] )

                T_SASI[r,m,mdot,rs] = tAd + tAc

                dataFileName \
                  = '.{:}_LegendrePowerSpectrum.dat'.format( ID[r,m,mdot,rs] )

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
for r in range( R.shape[0] ):
    for m in range( M.shape[0] ):
        for rs in range( Rs.shape[0] ):

            plotFileDirectory \
              = '/lump/data/accretionShockStudy/newRuns/newProductionRuns/{:}/'.format \
                ( ID[r,m,mdot,rs] )

            #if not isdir( plotFileDirectory ):
            #    print( '{:} does not exist. Skipping.' \
            #           .format( plotFileDirectory ) )
            #    continue

            tt  = t [r,m,mdot,rs]
            P1t = P1[r,m,mdot,rs]
            t0t = t0[r,m,mdot,rs]
            t1t = t1[r,m,mdot,rs]

            print( R[r], Rs[rs] )

            ind = np.where( ( tt >= t0t ) & ( tt <= t1t ) )[0]

            tF = tt[ind]

            logFt   = LogF  [r,m,mdot,rs]
            omegaRt = omegaR[r,m,mdot,rs]
            omegaIt = omegaI[r,m,mdot,rs]
            deltat  = delta [r,m,mdot,rs]

            F = FittingFunction \
                 ( tF - tF[0], logFt, \
                   omegaRt, omegaIt, deltat )

            tau = 1.0#T_SASI[0,m,mdot,rs]

            ind = np.where( ( tt < 710.0 ) & ( tt >= 0.0 ) )[0]
            #ind = np.where( tt < 100.0 )[0]

            if r == 1:
                ax.plot \
                  ( tt[ind]/tau, P1t[ind], 'g--', \
                    label = 'Rs={:} km (GR)'.format( Rs[rs] ) )

            else:
                ax.plot \
                  ( tt[ind]/tau, P1t[ind], '-', \
                    label = 'Rs={:} km'.format( Rs[rs] ) )
ax.grid()
ax.set_yscale( 'log' )
#ax.set_xlim( 0, 12 )
#xticks = np.linspace( 0, 12, 13, dtype = np.int64 )
#xticklabels = [ str( i ) for i in xticks ]
#ax.set_xticks( xticks )
#ax.set_xticklabels( xticklabels )
##ax.axvline( 150/T_SASI[0,0,0,0], label = r'$t=150\,\mathrm{ms}$' )

ax.set_title( r'$\texttt{{NR2D_M{:}_Rpns040}}$'.format( M[0] ), \
              fontsize = 15 )

ax.legend()
#fig.supxlabel( r'$t/T_{\mathrm{SASI,NR}}$' )
fig.supxlabel( 'Time [ms]' )
fig.supylabel( r'$H_{1}$ [cgs]' )
plt.subplots_adjust( hspace = 0.0, wspace = 0.3 )

#plt.savefig( '/home/kkadoogan/fig.PowerInLegendreMode.png', dpi = 300 )
plt.show()

import os
os.system( 'rm -rf __pycache__ ' )
