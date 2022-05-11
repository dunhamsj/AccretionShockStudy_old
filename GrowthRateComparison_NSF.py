#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( './Publication.sty' )

Root = '/lump/data/AccretionShockStudy/'

M    = np.array( [ '2.0' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '120' ], str )

G_GR     = np.loadtxt( 'G_GR.dat' )
G_err_GR = np.loadtxt( 'G_err_GR.dat' )
G_NR     = np.loadtxt( 'G_NR.dat' )
G_err_NR = np.loadtxt( 'G_err_NR.dat' )
T_GR     = np.loadtxt( 'T_GR.dat' )
T_err_GR = np.loadtxt( 'T_err_GR.dat' )
T_NR     = np.loadtxt( 'T_NR.dat' )
T_err_NR = np.loadtxt( 'T_err_NR.dat' )

M    = np.float64( M  )
Mdot = np.float64( Mdot )
Rs   = np.int64( Rs )

fig, ax = plt.subplots( 1, 1, figsize = (8,6) )

size = 15

ID_GR = 'GR2D_M2.0_Mdot0.3_Rs120'
Time, RsAve, RsMin, RsMax, P0, P1G, P2, P3, P4 \
  = np.loadtxt( '.' + ID_GR + '_DivV2_0.80-0.90_PowersInLegendreModes.dat' )
ID_NR = 'NR2D_M2.0_Mdot0.3_Rs120'
Time, RsAve, RsMin, RsMax, P0, P1N, P2, P3, P4 \
  = np.loadtxt( '.' + ID_NR + '_DivV2_0.80-0.90_PowersInLegendreModes.dat' )

ax.set_title( r'$M={:.1f}\ \mathrm{{M}}_{{\odot}},\ R_s={:d}\,\mathrm{{km}}$'.format( 2.0, 120 ), \
fontsize = size )
ax.plot( Time, P1N, c = 'r', label = 'NR' )
ax.plot( Time, P1G, c = 'b', label = 'GR' )
ax.text( 0.3, 0.9, r'$\tau_\mathrm{{NR}}: {:.3f}\ \mathrm{{ms}}$'.format \
             ( G_NR ), fontsize = size, transform = ax.transAxes, \
             color = 'red' )
ax.text( 0.3, 0.8, r'$\tau_\mathrm{{GR}}: {:.3f}\ \mathrm{{ms}}$'.format \
             ( G_GR ), fontsize = size, transform = ax.transAxes, \
             color = 'blue' )
ax.set_xlabel( 'Time [ms]', fontsize = size )
ax.set_yscale( 'log' )
ax.tick_params( axis = 'x', labelsize = size )
ax.tick_params( axis = 'y', labelsize = size )
ax.set_ylabel( 'Power [cgs]', fontsize = size )
ax.grid()
ax.legend(prop={'size':15})

#plt.show()
plt.savefig( 'fig.GrowthRateComparison.png', dpi = 300, \
             bbox_inches = 'tight' )

import os
os.system( 'rm -rf __pycache__ ' )
