#!/usr/bin/env python3

import numpy as np
from os.path import isdir
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from FitPowerToModel import FittingFunction

R    = 'NR'
M    = '2.8'
Mdot = np.array( [ '0.3', '0.5' ], str )
Rs   = '6.00e1'
#nX2  = np.array( [ '064', '128', '256' ], str )
nX2  = 'RPNS2.00e1'

arrShape = (Mdot.shape[0])

ID = np.empty( arrShape, object )

t  = np.empty( arrShape, object )
P0 = np.empty( arrShape, object )
P1 = np.empty( arrShape, object )
P2 = np.empty( arrShape, object )
P3 = np.empty( arrShape, object )
P4 = np.empty( arrShape, object )

fig, axs = plt.subplots( 1, 1 )

for N in range( Mdot.shape[0] ):

    ID[N] \
            = '{:}2D_M{:}_Mdot{:}_Rs{:}_RPNS2.00e1'.format \
        ( R, M, Mdot[N], Rs )

    print( ID[N] )
    plotFileDirectory \
      = '/lump/data/accretionShockStudy/{:}/'.format \
        ( ID[N] )

    if not isdir( plotFileDirectory ): continue

    plotFileBaseName = '{:}.plt'.format( ID[N] )

    dataFileName \
            = '.{:}_LegendrePowerSpectrum.dat'.format( ID[N] )

    t [N], P0[N], P1[N], P2[N], P3[N], P4[N] \
      = np.loadtxt( dataFileName )

    t_NR  = t [N]
    P0_NR = P0[N]
    P1_NR = P1[N]
    P2_NR = P2[N]
    P3_NR = P3[N]
    P4_NR = P4[N]

    ind = 2
#    axs.plot( t_NR[ind:], P0_NR[ind:], label = 'P0' )
    axs.plot( t_NR[ind:], P1_NR[ind:], label = 'Mdot = {:}'.format( Mdot[N] ) )
#    axs.plot( t_NR[ind:], P2_NR[ind:], label = 'P2' )
#    axs.plot( t_NR[ind:], P3_NR[ind:], label = 'P3' )
#    axs.plot( t_NR[ind:], P4_NR[ind:], label = 'P4' )

#    axs[N].set_xlim( 0.0, 100.0 )
    axs.set_yscale( 'log' )
    axs.set_title( 'M2.8_Rs6.00e1_RPNS2.00e1' )

#    if( N < Mdot.shape[0]-1 ): axs.xaxis.set_visible( False )

axs.legend()
fig.supxlabel( 'Time [ms]' )
fig.supylabel( r'$H_{1}$ [cgs]' )
plt.subplots_adjust( hspace = 0.0 )

plt.savefig( '/home/kkadoogan/fig.MdotComparison.png', dpi = 300 )
#plt.show()

import os
os.system( 'rm -rf __pycache__ ' )
