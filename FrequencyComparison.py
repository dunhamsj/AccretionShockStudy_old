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

T_SASI = np.empty( (M.shape[0],Rs.shape[0]), np.float64 )
T = np.array( [[25.07716059, 35.15626653, 55.64370891], \
               [19.88319758, 29.90565005, 50.06614409], \
               [19.91161272, 29.94128336, 41.91074087]], np.float64 )


for m in range( M.shape[0] ):
    for mdot in range( Mdot.shape[0] ):
        for rs in range( Rs.shape[0] ):

            ID \
              = 'GR1D_M{:}_Mdot{:}_Rs{:}'.format( M[m], Mdot[mdot], Rs[rs] )
            DataDirectory = Root + '{:}/'.format( ID )
            DataDirectory += '{:}.plt_00000000/'.format( ID )

            rInner = 4.0e1
            rOuter = np.float64( Rs[rs] )

            T_SASI[m,rs] \
              = TS.ComputeTimeScales( DataDirectory, rInner, rOuter )

M  = np.float64( M  )
Rs = np.int64( Rs )

## 1D, mass on x-axis
#
#c = [ 'r', 'b', 'm' ]
#s = [ '.', 's' ]
#
#for m in range( M.shape[0] ):
#    for rs in range( Rs.shape[0] ):
#        if m == 0:
#            plt.plot( M[m], T_SASI[m,rs], c[rs] + s[0], \
#                      label = r'$R_s={:d}$ km, est.'.format( Rs[rs] ) )
#            plt.plot( M[m], T[m,rs], c[rs] + s[1], \
#                      label = r'$R_s={:d}$ km, meas.'.format( Rs[rs] ) )
#        else:
#            plt.plot( M[m], T_SASI[m,rs], c[rs] + s[0] )
#            plt.plot( M[m], T[m,rs], c[rs] + s[1] )
#plt.xticks( M )
#plt.xlabel( r'$M_{\mathrm{PNS}}\ \left[\mathrm{M}_{\odot}\right]$' )
#plt.ylabel( r'$T_{\mathrm{SASI}}\ \left[\mathrm{ms}\right]$' )
#plt.legend()
#plt.show()
#plt.close()
#
## 1D, Rs on x-axis
#
#for m in range( M.shape[0] ):
#    for rs in range( Rs.shape[0] ):
#        if rs == 0:
#            plt.plot( Rs[rs], T_SASI[m,rs], c[m] + s[0], \
#                      label = r'$M={:.1f}\ M_\odot$, est.'.format( M[m] ) )
#            plt.plot( Rs[rs], T[m,rs], c[m] + s[1], \
#                      label = r'$M={:.1f}\ M_\odot$, meas.'.format( M[m] ) )
#        else:
#            plt.plot( Rs[rs], T_SASI[m,rs], c[m] + s[0] )
#            plt.plot( Rs[rs], T[m,rs], c[m] + s[1] )
#plt.xticks( Rs )
#plt.xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$' )
#plt.ylabel( r'$T_{\mathrm{SASI}}\ \left[\mathrm{ms}\right]$' )
#plt.legend()
#plt.show()
#plt.close()

# Heatmap

extent = [ 0, 3, 0, 3 ]

fig, ax = plt.subplots( 1, 1 )
fig.suptitle( 'GR' )

im = ax.imshow( ( T_SASI - T ) / T_SASI, \
                origin = 'lower', \
                extent = extent, \
                cmap = 'RdBu', \
                interpolation = 'bilinear', \
                aspect = 'auto' )
yticklabels = M
xticklabels = Rs
xticks = np.array( [ 1.0, 3.0, 5.0 ], np.float64 ) / 2.0
yticks = np.array( [ 1.0, 3.0, 5.0 ], np.float64 ) / 2.0
ax.set_xticks( xticks )
ax.set_yticks( yticks )
ax.set_xticklabels( xticklabels )
ax.set_yticklabels( yticklabels )

ax.set_xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$' )
ax.set_ylabel( r'$M\ \left[\mathrm{M}_{\odot}\right]$' )

cbar = fig.colorbar( im )
cbar.set_label( r'$( T_{\mathrm{SASI,est}} - T_{\mathrm{SASI,meas}} ) / T_{\mathrm{SASI,est}}$' )

plt.show()
#plt.savefig( 'fig.FrequencyComparison.png', dpi = 300 )
