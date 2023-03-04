#!/usr/bin/env python3

import numpy as np
from os.path import isdir
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from FitPowerToModel import FittingFunction

R    = np.array( [ 'NR', 'GR' ], str )
M    = np.array( [ '2.8', '2.8' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '6.00e1', '7.50e1', '9.00e1' ], str )

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

for r in range( R.shape[0] ):
    for m in range( M.shape[0] ):
        for mdot in range( Mdot.shape[0] ):
            for rs in range( Rs.shape[0] ):

                ID[r,m,mdot,rs] \
                  = '{:}2D_M{:}_Mdot{:}_Rs{:}_RPNS2.00e1'.format \
                     ( R[r], M[m], Mdot[mdot], Rs[rs] )

                plotFileDirectory \
                  = '/lump/data/accretionShockStudy/{:}/'.format \
                    ( ID[r,m,mdot,rs] )

                if not isdir( plotFileDirectory ): continue

                plotFileBaseName = '{:}.plt'.format( ID[r,m,mdot,rs] )

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

fig, axs = plt.subplots( Rs.shape[0], figsize = (12,9) )

nr   = 0
gr   = 1
mdot = 0
for m in range( M.shape[0] ):
    if m > 0: continue
    for rs in range( Rs.shape[0] ):

        t_NR = t[nr,m,mdot,rs]
        t_GR = t[gr,m,mdot,rs]

        P1_NR = P1[nr,m,mdot,rs]
        P1_GR = P1[gr,m,mdot,rs]

        t0_NR = t0[nr,m,mdot,rs]
        t0_GR = t0[gr,m,mdot,rs]

        t1_NR = t1[nr,m,mdot,rs]
        t1_GR = t1[gr,m,mdot,rs]

        ind_NR = np.where( ( t_NR >= t0_NR ) & ( t_NR <= t1_NR ) )[0]
        ind_GR = np.where( ( t_GR >= t0_GR ) & ( t_GR <= t1_GR ) )[0]

        tF_NR = t_NR[ind_NR]
        tF_GR = t_GR[ind_GR]

        logF_NR = LogF[nr,m,mdot,rs]
        logF_GR = LogF[gr,m,mdot,rs]

        omegaR_NR = omegaR[nr,m,mdot,rs]
        omegaR_GR = omegaR[gr,m,mdot,rs]

        omegaI_NR = omegaI[nr,m,mdot,rs]
        omegaI_GR = omegaI[gr,m,mdot,rs]

        delta_NR = delta[nr,m,mdot,rs]
        delta_GR = delta[gr,m,mdot,rs]

        F_NR = FittingFunction \
                 ( tF_NR - tF_NR[0], logF_NR, omegaR_NR, omegaI_NR, delta_NR )
        F_GR = FittingFunction \
                 ( tF_GR - tF_GR[0], logF_GR, omegaR_GR, omegaI_GR, delta_GR )

        axs[rs].plot( t_NR , P1_NR )
        axs[rs].plot( t_GR , P1_GR )
        axs[rs].plot( tF_NR, np.exp( F_NR ) )
        axs[rs].plot( tF_GR, np.exp( F_GR ) )

        axs[rs].text( 0.1, 0.7, \
        r'$\omega^{{\mathrm{{NR}}}}={:.3f}\ \mathrm{{ms}}^{{-1}}$' \
                      .format( omegaR_NR), \
                      transform = axs[rs].transAxes )
        axs[rs].text( 0.1, 0.8, \
        r'$\omega^{{\mathrm{{GR}}}}={:.3f}\ \mathrm{{ms}}^{{-1}}$' \
                      .format( omegaR_GR), \
                      transform = axs[rs].transAxes )

        axs[rs].text( 0.23, 0.7, \
        r'$T^{{\mathrm{{NR}}}}={:.3f}\ \mathrm{{ms}}$' \
                      .format( 2.0 * np.pi / omegaI_NR), \
                      transform = axs[rs].transAxes )
        axs[rs].text( 0.23, 0.8, \
        r'$T^{{\mathrm{{GR}}}}={:.3f}\ \mathrm{{ms}}$' \
                      .format( 2.0 * np.pi / omegaI_GR), \
                      transform = axs[rs].transAxes )

        axs[rs].text( 0.1, 0.9, 'M{:}_Rs{:}'.format( M[m], Rs[rs] ), \
                      transform = axs[rs].transAxes )

        axs[rs].set_xlim( 0.0, 100.0 )
        axs[rs].set_yscale( 'log' )

        if( rs < Rs.shape[0]-1 ): axs[rs].xaxis.set_visible( False )

fig.supxlabel( 'Time [ms]' )
fig.supylabel( r'$H_{1}$ [cgs]' )
plt.subplots_adjust( hspace = 0.0, wspace = 0.3 )

#plt.savefig( '/home/kkadoogan/fft_comp_LateStage.png', dpi = 300 )
plt.show()

import os
os.system( 'rm -rf __pycache__ ' )
