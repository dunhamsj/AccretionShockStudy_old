#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( './Publication.sty' )

Root = '/lump/data/AccretionShockStudy/'

M    = np.array( [ '1.4', '2.0', '2.8' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '120', '150', '180' ], str )

Field = 'DivV2'

G_GR     = np.loadtxt( 'G_GR_{:}.dat'    .format( Field ) )
G_err_GR = np.loadtxt( 'G_err_GR_{:}.dat'.format( Field ) )
G_NR     = np.loadtxt( 'G_NR_{:}.dat'    .format( Field ) )
G_err_NR = np.loadtxt( 'G_err_NR_{:}.dat'.format( Field ) )

M    = np.float64( M    )
Mdot = np.float64( Mdot )
Rs   = np.int64  ( Rs   )

# Heatmap

fig, axs = plt.subplots( 2, 1, figsize = (9,12) )

size = 15

ID_GR = 'GR2D_M2.8_Mdot0.3_Rs180'
Time, RsAve, RsMin, RsMax, P0, P1G, P2, P3, P4 \
  = np.loadtxt( '.' + ID_GR + '_DivV2_0.80-0.90_PowersInLegendreModes.dat' )
ID_NR = 'NR2D_M2.8_Mdot0.3_Rs180'
Time, RsAve, RsMin, RsMax, P0, P1N, P2, P3, P4 \
  = np.loadtxt( '.' + ID_NR + '_DivV2_0.80-0.90_PowersInLegendreModes.dat' )

ind = np.where( Time < 150.0 )[0]

axs[0].set_title( 'M2.8_Rs180', fontsize = size )
axs[0].plot( Time[ind], P1N[ind], c = 'r', label = 'NR' )
axs[0].plot( Time[ind], P1G[ind], c = 'b', label = 'GR' )
axs[0].text( 0.3, 0.9, r'$\omega_\mathrm{{NR}}: {:.3f}\ \mathrm{{Hz}}$'.format \
             ( G_NR[2,2] ), fontsize = size, transform = axs[0].transAxes, \
             color = 'red' )
axs[0].text( 0.3, 0.8, r'$\omega_\mathrm{{GR}}: {:.3f}\ \mathrm{{Hz}}$'.format \
             ( G_GR[2,2] ), fontsize = size, transform = axs[0].transAxes, \
             color = 'blue' )
axs[0].set_xlabel( 'Time [ms]', fontsize = size )
axs[0].set_yscale( 'log' )
axs[0].tick_params( axis = 'x', labelsize = size )
axs[0].tick_params( axis = 'y', labelsize = size )
axs[0].set_ylabel( 'Power [cgs]', fontsize = size )
axs[0].grid()
axs[0].legend(prop={'size':15})

extent = [ 0, Rs.shape[0], 0, M.shape[0] ]

Data = ( G_GR - G_NR ) / G_GR
im = axs[1].imshow( Data, \
                    origin = 'lower', \
                    extent = extent, \
                    cmap = 'viridis', \
#                    vmin = -0.01, vmax = +0.01, \
                    interpolation = 'none', \
                    aspect = 'auto' )

xticks = np.empty( Rs.shape[0], np.float64 )
yticks = np.empty( M.shape [0], np.float64 )
for i in range( xticks.shape[0] ):
    j = i + 1
    xticks[i] = ( 2.0 * np.float64( j ) - 1.0 ) / 2.0
for i in range( yticks.shape[0] ):
    j = i + 1
    yticks[i] = ( 2.0 * np.float64( j ) - 1.0 ) / 2.0

for i in range( xticks.shape[0] ):
    for j in range( yticks.shape[0] ):
        axs[1].text( xticks[i]-0.15, yticks[j], \
                     '{:.3e}'.format( Data[j,i] ), c = 'w' )

axs[1].set_xticks( xticks )
axs[1].set_yticks( yticks )
axs[1].set_xticklabels( Rs, size = size )
axs[1].set_yticklabels( M , size = size )

axs[1].set_xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$'   , fontsize = size )
axs[1].set_ylabel( r'$M\ \left[\mathrm{M}_{\odot}\right]$', fontsize = size )

cbar = fig.colorbar( im, pad = 0.01, location='top' )
cbar.set_label( r'$( \omega_{\mathrm{GR}} - \omega_{\mathrm{NR}} ) / \omega_{\mathrm{GR}}$', \
                size = size+2, labelpad = 10 )
cbar.ax.tick_params( labelsize = size )
plt.subplots_adjust( hspace = 0.3 )
plt.show()
#plt.savefig( \
#  '/home/kkadoogan/fig.GrowthRateComparison_HeatMap_{:}.png'.format( Field ), \
#  dpi = 300, bbox_inches = 'tight' )

import os
os.system( 'rm -rf __pycache__ ' )
