#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( './Publication.sty' )

from TimeScales import TimeScales

Root = '/lump/data/accretionShockStudy/'

M    = np.array( [ '1.4', '2.0', '2.8' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '120', '150', '180' ], str )

TS = TimeScales()

T_Mul    = np.empty( (M.shape[0],Rs.shape[0]), np.float64 )
T_GR     = np.loadtxt( 'T_GR_DivV2.dat' )
T_err_GR = np.loadtxt( 'T_err_GR_DivV2.dat' )
T_NR     = np.loadtxt( 'T_NR_DivV2.dat' )
T_err_NR = np.loadtxt( 'T_err_NR_DivV2.dat' )

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
            ax.plot    ( Rs[rs], T_Mul  [m,rs], c[m] + s[0], \
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
            ax.plot    ( Rs[rs], T_Mul  [m,rs], c[m] + s[0] )
            ax.errorbar( Rs[rs], T_GR   [m,rs], \
                         yerr = T_err_GR[m,rs], \
                         fmt = c[m] + s[1], markerfacecolor='none' )
            ax.errorbar( Rs[rs], T_NR   [m,rs], \
                         yerr = T_err_NR[m,rs], \
                         fmt = c[m] + s[2], markerfacecolor='none' )
handles, labels = ax.get_legend_handles_labels()

h = handles
l = labels
h[0], l[0] = handles[0], labels[0]
h[1], l[1] = handles[3], labels[3]
h[2], l[2] = handles[6], labels[6]
h[3], l[3] = handles[1], labels[1]
h[4], l[4] = handles[4], labels[4]
h[5], l[5] = handles[7], labels[7]
h[6], l[6] = handles[2], labels[2]
h[7], l[7] = handles[5], labels[5]
h[8], l[8] = handles[8], labels[8]
ax.legend(h,l)

ax.set_xticks( Rs )
ax.set_xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$' )
ax.set_ylabel( r'$T_{\mathrm{SASI}}\ \left[\mathrm{ms}\right]$' )
#plt.show()
plt.savefig( '/home/kkadoogan/fig.FrequencyComparison_1D_Rs.png', dpi = 300 )
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

# Heatmap

extent = [ 0, Rs.shape[0], 0, M.shape[0] ]

fig, ax = plt.subplots( 1, 1 )

Data = ( T_NR - T_GR ) / T_NR

print( Data.max() )
print( Data.min() )
im = ax.imshow( Data, \
                origin = 'lower', \
                extent = extent, \
                cmap = 'RdBu', \
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
ax.set_xticks( xticks )
ax.set_yticks( yticks )
ax.set_xticklabels( xticklabels )
ax.set_yticklabels( yticklabels )

ax.set_xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$' )
ax.set_ylabel( r'$M\ \left[\mathrm{M}_{\odot}\right]$' )

cbar = fig.colorbar( im )
cbar.set_label( r'$( T_{\mathrm{GR}} - T_{\mathrm{NR}} ) / T_{\mathrm{GR}}$' )
#plt.show()
plt.savefig( '/home/kkadoogan/fig.FrequencyComparison_HeatMap.png', dpi = 300 )

import os
os.system( 'rm -rf __pycache__ ' )
