#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
plt.style.use( './Publication.sty' )

from TimeScales import TimeScales

Root = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Euler_Relativistic_IDEAL/'

M    = np.linspace( 1.0, 3.0, 11 )
Mdot = np.array( [ 0.3 ], np.float64 )
Rs   = np.linspace( 110, 200, 10 )
T_GR = np.loadtxt( 'T_GR.dat' )

TS = TimeScales()

T_Mul = np.empty( (M.shape[0],Rs.shape[0]), np.float64 )

#for m in range( M.shape[0] ):
#    for mdot in range( Mdot.shape[0] ):
#        for rs in range( Rs.shape[0] ):
#
#            ID \
#              = 'GR1D_M{:.1f}_Mdot{:.1f}_Rs{:g}'.format \
#                  ( M[m], Mdot[mdot], Rs[rs] )
#            DataDirectory = Root + '{:}.plt_00000000/'.format( ID )
#
#            rInner = 4.0e1
#            rOuter = np.float64( Rs[rs] )
#
#            T_Mul[m,rs] = TS.ComputeTimeScales( DataDirectory, rInner, rOuter )
#
#np.savetxt( 'T_Mul.dat', T_Mul )
T_Mul = np.loadtxt( 'T_Mul.dat' )
size = 10
fig, ax = plt.subplots( 1, 1 )
extent = [ Rs.min()-5, Rs.max()+5, M.min()-0.1, M.max()+0.1 ]
vmin = 10
vmax = 80
im = ax.imshow( T_Mul, \
                origin = 'lower', \
                extent = extent, \
                cmap = 'viridis', \
                vmin = vmin, vmax = vmax, \
                interpolation = 'bilinear', \
                aspect = 'auto' )

# from: https://matplotlib.org/3.5.0/gallery/images_contours_and_fields/
#       contour_label_demo.html
def fmt(x):
    s = f"{x:g} ms"
    return s

contour = ax.contour( T_Mul, extent = extent, colors = 'white' )
ax.clabel( contour, contour.levels, inline=True, fmt=fmt, fontsize=10)

xticklabels = Rs
yticklabels = M
xticks = Rs
yticks = M

ax.set_xticks( xticks )
ax.set_yticks( yticks )
ax.set_xticklabels( xticklabels, size = size )
ax.set_yticklabels( yticklabels, size = size )
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax.set_xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$'   , fontsize = size )
ax.set_ylabel( r'$M\ \left[\mathrm{M}_{\odot}\right]$', fontsize = size )

cbar = fig.colorbar( im, pad = 0.01, location='top' )
cbar.set_label( r'$T_{\mathrm{SASI}}$ [ms]', size = size+2, labelpad = 10 )
cbar.ax.tick_params( labelsize = size )

M = np.array( [ 1.4, 2.0, 2.8 ], np.float64 )
Rs = np.array( [ 120.0, 150.0, 180.0 ], np.float64 )
for m in range( M.shape[0] ):
    for rs in range( Rs.shape[0] ):
        ax.scatter( Rs[rs], M[m], c=T_GR[m,rs], vmin = vmin, vmax = vmax, \
                    edgecolors = 'k' )

#plt.show()
plt.savefig( 'fig.T_SASI.png', dpi = 300 )
import os
os.system( 'rm -rf __pycache__ ' )
