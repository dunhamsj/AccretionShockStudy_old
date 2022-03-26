#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( './Publication.sty' )

from TimeScales import TimeScales

Root = '/lump/data/AccretionShockStudy/'

M    = np.array( [ '1.4', '2.0', '2.8' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '120', '150', '180' ], str )

TS = TimeScales()

T_Mul    = np.empty( (M.shape[0],Rs.shape[0]), np.float64 )
T_GR     = np.loadtxt( 'T_GR.dat' )
T_err_GR = np.loadtxt( 'T_err_GR.dat' )
T_NR     = np.loadtxt( 'T_NR.dat' )
T_err_NR = np.loadtxt( 'T_err_NR.dat' )

for m in range( M.shape[0] ):
    for mdot in range( Mdot.shape[0] ):
        for rs in range( Rs.shape[0] ):

            ID \
              = 'GR1D_M{:}_Mdot{:}_Rs{:}'.format( M[m], Mdot[mdot], Rs[rs] )
            DataDirectory = Root + '{:}/'.format( ID )
            DataDirectory += '{:}.plt_00000000/'.format( ID )

            rInner = 4.0e1
            rOuter = np.float64( Rs[rs] )

            T_Mul[m,rs] = TS.ComputeTimeScales( DataDirectory, rInner, rOuter )

M  = np.float64( M  )
Rs = np.int64( Rs )

# 1D
c = [ 'r', 'b', 'm' ]
s = [ 'x', '.', 's' ]

# Rs on x-axis

NX = 5.0
NY = NX * 3/4
fig, ax = plt.subplots( 1, 1, figsize = (NX,NY) )

for rs in range( Rs.shape[0] ):
    for m in range( M.shape[0] ):
        if rs == 0:
            ax.plot( Rs[rs], T_Mul[m,rs], c[m] + s[0], \
                     label = r'$M={:.1f}\ M_\odot$, Mul.'.format( M[m] ) )
            ax.errorbar( Rs[rs], T_GR   [m,rs], \
                         yerr = T_err_GR[m,rs], \
                         fmt = c[m] + s[1], markerfacecolor='none', \
                         label = r'$M={:.1f}\ M_\odot$, GR'.format( M[m] ) )
            ax.errorbar( Rs[rs], T_NR   [m,rs], \
                         yerr = T_err_NR[m,rs], \
                         fmt = c[m] + s[2], markerfacecolor='none' , \
                         label = r'$M={:.1f}\ M_\odot$, NR'.format( M[m] ) )
        else:
            ax.plot( Rs[rs], T_Mul[m,rs], c[m] + s[0] )
            ax.errorbar( Rs[rs], T_GR   [m,rs], \
                         yerr = T_err_GR[m,rs], \
                         fmt = c[m] + s[1], markerfacecolor='none' )
            ax.errorbar( Rs[rs], T_NR   [m,rs], \
                         yerr = T_err_NR[m,rs], \
                         fmt = c[m] + s[2], markerfacecolor='none' )
ax.set_xticks( Rs )
ax.set_xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$' )
ax.set_ylabel( r'$T_{\mathrm{SASI}}\ \left[\mathrm{ms}\right]$' )
ax.legend()
#plt.show()
plt.savefig( 'fig.FrequencyComparison_1D_Rs.png', dpi = 300 )
plt.close()

## mass on x-axis
#
#for m in range( M.shape[0] ):
#    for rs in range( Rs.shape[0] ):
#        if m == 0:
#            plt.plot( M[m], T_Mul[m,rs], c[rs] + s[0], \
#                      label = r'$R_s={:d}$ km, Mul.'.format( Rs[rs] ) )
#            plt.errorbar( M[m], T_GR     [m,rs], \
#                          yerr = T_err_GR[m,rs], \
#                          fmt = c[rs] + s[1], \
#                      label = r'$R_s={:d}$ km, meas.'.format( Rs[rs] ) )
#        else:
#            plt.plot( M[m], T_Mul[m,rs], c[rs] + s[0] )
#            plt.errorbar( M[m], T_GR     [m,rs], \
#                          yerr = T_err_GR[m,rs], \
#                          fmt = c[rs] + s[1] )
#plt.xticks( M )
#plt.xlabel( r'$M_{\mathrm{PNS}}\ \left[\mathrm{M}_{\odot}\right]$' )
#plt.ylabel( r'$T_{\mathrm{SASI}}\ \left[\mathrm{ms}\right]$' )
#plt.legend()
##plt.show()
#plt.savefig( 'fig.FrequencyComparison_1D_Mass.png', dpi = 300 )
#plt.close()

## Heatmap
#
#extent = [ 0, Rs.shape[0], 0, M.shape[0] ]
#
#fig, ax = plt.subplots( 1, 1 )
#
#im = ax.imshow( ( T_GR - T_GR ) / T_GR, \
#                origin = 'lower', \
#                extent = extent, \
#                cmap = 'RdBu', \
#                vmin = -0.01, vmax = +0.01, \
#                interpolation = 'bilinear', \
#                aspect = 'auto' )
#yticklabels = M
#xticklabels = Rs
#xticks = np.empty( Rs.shape[0], np.float64 )
#yticks = np.empty( M.shape [0], np.float64 )
#for i in range( xticks.shape[0] ):
#    j = i + 1
#    xticks[i] = ( 2.0 * np.float64( j ) - 1.0 ) / 2.0
#for i in range( yticks.shape[0] ):
#    j = i + 1
#    yticks[i] = ( 2.0 * np.float64( j ) - 1.0 ) / 2.0
#ax.set_xticks( xticks )
#ax.set_yticks( yticks )
#ax.set_xticklabels( xticklabels )
#ax.set_yticklabels( yticklabels )
#
#ax.set_xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$' )
#ax.set_ylabel( r'$M\ \left[\mathrm{M}_{\odot}\right]$' )
#
#cbar = fig.colorbar( im )
#cbar.set_label( r'$( T_{\mathrm{GR}} - T_{\mathrm{NR}} ) / T_{\mathrm{GR}}$' )
##plt.show()
#plt.savefig( 'fig.FrequencyComparison_HeatMap.png', dpi = 300 )

import os
os.system( 'rm -rf __pycache__ ' )
