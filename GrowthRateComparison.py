#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( './Publication.sty' )

Root = '/lump/data/AccretionShockStudy/'

M    = np.array( [ '1.4', '2.0' ], str )
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

# 1D
c = [ 'r', 'b', 'm' ]
s = [ 'x', '.', 's' ]

# Rs on x-axis

for m in range( M.shape[0] ):
    for rs in range( Rs.shape[0] ):
        if rs == 0:
            plt.errorbar( Rs[rs], G_GR   [m,rs], \
                          yerr = G_err_GR[m,rs], \
                          fmt = c[m] + s[1], markerfacecolor='none', \
                      label = r'$M={:.1f}\ M_\odot$, GR'.format( M[m] ) )
            plt.errorbar( Rs[rs], G_NR   [m,rs], \
                          yerr = G_err_NR[m,rs], \
                          fmt = c[m] + s[2], markerfacecolor='none' , \
                      label = r'$M={:.1f}\ M_\odot$, NR'.format( M[m] ) )
        else:
            plt.errorbar( Rs[rs], G_GR   [m,rs], \
                          yerr = G_err_GR[m,rs], \
                          fmt = c[m] + s[1], markerfacecolor='none' )
            plt.errorbar( Rs[rs], G_NR   [m,rs], \
                          yerr = G_err_NR[m,rs], \
                          fmt = c[m] + s[2], markerfacecolor='none' )
plt.xticks( Rs )
plt.xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$' )
plt.ylabel( r'$\tau_{\mathrm{SASI}}\ \left[\mathrm{ms}\right]$' )
plt.legend()
#plt.show()
plt.savefig( 'fig.GrowthRateComparison_1D_Rs.png', dpi = 300 )
plt.close()

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

fig, ax = plt.subplots( 1, 1 )
#fig.suptitle( 'GR vs. NR' )

im = ax.imshow( ( G_GR - G_NR ) / G_GR, \
                origin = 'lower', \
                extent = extent, \
                cmap = 'viridis', \
#                vmin = -0.01, vmax = +0.01, \
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

size = 15

ax.set_xticks( xticks )
ax.set_yticks( yticks )
ax.set_xticklabels( xticklabels, size = size )
ax.set_yticklabels( yticklabels, size = size )

ax.set_xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$'   , fontsize = size )
ax.set_ylabel( r'$M\ \left[\mathrm{M}_{\odot}\right]$', fontsize = size )

cbar = fig.colorbar( im )
cbar.set_label( r'$( \tau_{\mathrm{GR}} - \tau_{\mathrm{NR}} ) / \tau_{\mathrm{GR}}$' )
#plt.show()
plt.savefig( 'fig.GrowthRateComparison_HeatMap.png', dpi = 300 )

import os
os.system( 'rm -rf __pycache__ ' )
