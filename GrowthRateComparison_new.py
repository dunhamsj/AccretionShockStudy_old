#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( './Publication.sty' )

Root = '/lump/data/accretionShockStudy/'

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

fig, ax = plt.subplots( 1, 1, figsize = (9,12) )

extent = [ 0, Rs.shape[0], 0, M.shape[0] ]

cmap = plt.get_cmap( 'viridis' )

#Data = ( G_GR - G_NR ) / G_GR
Data = G_GR
im = ax.imshow( Data, \
                    origin = 'lower', \
                    extent = extent, \
                    cmap = cmap, \
#                    vmin = -0.01, vmax = +0.01, \
                    interpolation = 'bilinear', \
                    aspect = 'auto' )

a1 = min( Data.min(), 26.687800769228687 )
a2 = max( Data.max(), 31.128138422245645 )
norm = plt.Normalize( a1, a2 )

c1 = cmap( norm( 31.128138422245645 ) )
ax.plot( 2.0, 2.5, 'o', color = c1, markeredgecolor = 'k' )

c2 = cmap( norm( 26.687800769228687 ) )
ax.plot( 2.5, 2.0, 'o', color = c2, markeredgecolor = 'k' )

c3 = cmap( norm( 27.2280909909708 ) )
ax.plot( 2.0, 2.0, 'o', color = c3, markeredgecolor = 'k' )

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
        ax.text( xticks[i]-0.15, yticks[j], \
                     '{:.3e}'.format( Data[j,i] ), c = 'w' )

size = 15

ax.set_xticks( xticks )
ax.set_yticks( yticks )
ax.set_xticklabels( Rs, size = size )
ax.set_yticklabels( M , size = size )

ax.set_xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$'   , fontsize = size )
ax.set_ylabel( r'$M\ \left[\mathrm{M}_{\odot}\right]$', fontsize = size )

cbar = fig.colorbar( im, pad = 0.01, location='top' )
cbar.set_label( r'$\omega_{\mathrm{GR}}\ \left[\mathrm{MHz}\right]$', \
                size = size+2, labelpad = 10 )
cbar.ax.tick_params( labelsize = size )
#plt.show()
plt.savefig( \
  '/home/kkadoogan/fig.GrowthRateComparison_HeatMap_{:}.png'.format( Field ), \
  dpi = 300, bbox_inches = 'tight' )

import os
os.system( 'rm -rf __pycache__ ' )
