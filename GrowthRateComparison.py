#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( './Publication.sty' )

Root = '/lump/data/AccretionShockStudy/'

M    = np.array( [ '1.4', '2.0', '2.8' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '120', '150', '180' ], str )

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

## 1D
#c = [ 'r', 'b', 'm' ]
#s = [ 'x', '.', 's' ]
#
## Rs on x-axis
#
#for m in range( M.shape[0] ):
#    for rs in range( Rs.shape[0] ):
#        if rs == 0:
#            plt.errorbar( Rs[rs], G_GR   [m,rs], \
#                          yerr = G_err_GR[m,rs], \
#                          fmt = c[m] + s[1], markerfacecolor='none', \
#                      label = r'$M={:.1f}\ M_\odot$, GR'.format( M[m] ) )
#            plt.errorbar( Rs[rs], G_NR   [m,rs], \
#                          yerr = G_err_NR[m,rs], \
#                          fmt = c[m] + s[2], markerfacecolor='none' , \
#                      label = r'$M={:.1f}\ M_\odot$, NR'.format( M[m] ) )
#        else:
#            plt.errorbar( Rs[rs], G_GR   [m,rs], \
#                          yerr = G_err_GR[m,rs], \
#                          fmt = c[m] + s[1], markerfacecolor='none' )
#            plt.errorbar( Rs[rs], G_NR   [m,rs], \
#                          yerr = G_err_NR[m,rs], \
#                          fmt = c[m] + s[2], markerfacecolor='none' )
#plt.xticks( Rs )
#plt.xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$' )
#plt.ylabel( r'$\tau_{\mathrm{SASI}}\ \left[\mathrm{ms}\right]$' )
#plt.legend()
##plt.show()
#plt.savefig( 'fig.GrowthRateComparison_1D_Rs.png', dpi = 300 )
#plt.close()

## mass on x-axis
#
#for m in range( M.shape[0] ):
#    for rs in range( Rs.shape[0] ):
#        if m == 0:
#            plt.errorbar( M[m], G_GR     [m,rs], \
#                          yerr = G_err_GR[m,rs], \
#                          fmt = c[rs] + s[1], \
#                      label = r'$R_s={:d}$ km, meas.'.format( Rs[rs] ) )
#        else:
#            plt.errorbar( M[m], G_GR     [m,rs], \
#                          yerr = G_err_GR[m,rs], \
#                          fmt = c[rs] + s[1] )
#plt.xticks( M )
#plt.xlabel( r'$M_{\mathrm{PNS}}\ \left[\mathrm{M}_{\odot}\right]$' )
#plt.ylabel( r'$\tau_{\mathrm{SASI}}\ \left[\mathrm{ms}\right]$' )
#plt.legend()
##plt.show()
#plt.savefig( 'fig.GrowthRateComparison_1D_Mass.png', dpi = 300 )
#plt.close()

# Heatmap

extent = [ 0, Rs.shape[0], 0, M.shape[0] ]

fig, axs = plt.subplots( 2, 1, figsize = (9,12) )

size = 15

ID_GR = 'GR2D_M2.0_Mdot0.3_Rs120'
Time, RsAve, RsMin, RsMax, P0, P1G, P2, P3, P4 \
  = np.loadtxt( '.' + ID_GR + '_DivV2_0.80-0.90_PowersInLegendreModes.dat' )
ID_NR = 'NR2D_M2.0_Mdot0.3_Rs120'
Time, RsAve, RsMin, RsMax, P0, P1N, P2, P3, P4 \
  = np.loadtxt( '.' + ID_NR + '_DivV2_0.80-0.90_PowersInLegendreModes.dat' )

axs[0].set_title( 'M2.0_Rs120', fontsize = size )
axs[0].plot( Time, P1N, c = 'r', label = 'NR' )
axs[0].plot( Time, P1G, c = 'b', label = 'GR' )
axs[0].text( 0.3, 0.9, r'$\tau_\mathrm{{NR}}: {:.3f}\ \mathrm{{ms}}$'.format \
             ( G_NR[1,0] ), fontsize = size, transform = axs[0].transAxes, \
             color = 'red' )
axs[0].text( 0.3, 0.8, r'$\tau_\mathrm{{GR}}: {:.3f}\ \mathrm{{ms}}$'.format \
             ( G_GR[1,0] ), fontsize = size, transform = axs[0].transAxes, \
             color = 'blue' )
axs[0].set_xlabel( 'Time [ms]', fontsize = size )
axs[0].set_yscale( 'log' )
axs[0].tick_params( axis = 'x', labelsize = size )
axs[0].tick_params( axis = 'y', labelsize = size )
axs[0].set_ylabel( 'Power [cgs]', fontsize = size )
axs[0].grid()
axs[0].legend(prop={'size':15})

im = axs[1].imshow( ( G_GR - G_NR ) / G_GR, \
                    origin = 'lower', \
                    extent = extent, \
                    cmap = 'viridis', \
#                    vmin = -0.01, vmax = +0.01, \
                    interpolation = 'bilinear', \
                    aspect = 'auto' )
yticklabels = M
xticklabels = Rs
xticks = np.empty( Rs.shape[0], np.float64 )
yticks = np.empty( M.shape [0], np.float64 )
for i in range( xticks.shape[0] ):
    j = i + 1
    xticks[i] = ( 2.0 * np.float64( j ) - 1.0 ) / 2.0
for i in range( yticks.shape[0] ):
    j = i + 1
    yticks[i] = ( 2.0 * np.float64( j ) - 1.0 ) / 2.0

axs[1].set_xticks( xticks )
axs[1].set_yticks( yticks )
axs[1].set_xticklabels( xticklabels, size = size )
axs[1].set_yticklabels( yticklabels, size = size )

axs[1].set_xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$'   , fontsize = size )
axs[1].set_ylabel( r'$M\ \left[\mathrm{M}_{\odot}\right]$', fontsize = size )

cbar = fig.colorbar( im, pad = 0.01, location='top' )
cbar.set_label( r'$( \tau_{\mathrm{GR}} - \tau_{\mathrm{NR}} ) / \tau_{\mathrm{GR}}$', size = size+2, labelpad = 10 )
cbar.ax.tick_params( labelsize = size )
plt.subplots_adjust( hspace = 0.3 )
#plt.show()
plt.savefig( 'fig.GrowthRateComparison_HeatMap.png', dpi = 300, \
             bbox_inches = 'tight' )

import os
os.system( 'rm -rf __pycache__ ' )
