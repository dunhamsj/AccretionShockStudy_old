#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( './Publication.sty' )

from TimeScales import TimeScales
from UtilitiesModule import GetNorm

rootDirectory = '/lump/data/accretionShockStudy/'

M    = np.array( [ '1.4', '2.0', '2.8' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '120', '150', '180' ], str )

Field = 'DivV2'

T_GR     = np.loadtxt( 'T_GR_{:}.dat'    .format( Field ) )
T_err_GR = np.loadtxt( 'T_err_GR_{:}.dat'.format( Field ) )
T_NR     = np.loadtxt( 'T_NR_{:}.dat'    .format( Field ) )
T_err_NR = np.loadtxt( 'T_err_NR_{:}.dat'.format( Field ) )

TS = TimeScales()
T_Mul = np.empty( (M.shape[0],Rs.shape[0]), np.float64 )
for m in range( M.shape[0] ):
    for mdot in range( Mdot.shape[0] ):
        for rs in range( Rs.shape[0] ):

            ID \
              = 'GR1D_M{:}_Mdot{:}_Rs{:}'.format( M[m], Mdot[mdot], Rs[rs] )
            DataDirectory = rootDirectory + '{:}/'.format( ID )
            DataDirectory += '{:}.plt_00000000/'.format( ID )

            rInner = 4.0e1
            rOuter = np.float64( Rs[rs] )

            T_Mul[m,rs] \
              = TS.ComputeTimeScales( DataDirectory, rInner, rOuter )

dT = np.max( np.abs( ( T_GR - T_NR ) / T_NR ) )
print( '|GR-NR|/|NR|: {:.3e}'.format( dT ) )
dT = np.max( np.abs( ( T_GR - T_Mul ) / T_Mul ) )
print( '|GR-Mul|/|Mul|: {:.3e}'.format( dT ) )
dT = np.max( np.abs( ( T_NR - T_Mul ) / T_Mul ) )
print( '|NR-Mul|/|Mul|: {:.3e}'.format( dT ) )

M    = np.float64( M    )
Mdot = np.float64( Mdot )
Rs   = np.int64  ( Rs   )

c = [ 'r', 'b', 'm' ]
s = [ 'x', '.', 's' ]

NX = 5.0
NY = NX * 3/4
fig, ax = plt.subplots( 1, 1, figsize = (NX,NY) )

for rs in range( Rs.shape[0] ):
    for m in range( M.shape[0] ):

        if rs == 0:

            ax.plot    ( Rs[rs], T_Mul  [m,rs], c[m] + s[0], \
                         label \
                         = r'M$\mathrm{{{{\"u}}}}$lEst_M{:.1f}_Rs{:d}'.format \
                           ( M[m], Rs[rs] ) )
            ax.errorbar( Rs[rs], T_GR   [m,rs], \
                         yerr = T_err_GR[m,rs], \
                         fmt = c[m] + s[1], markerfacecolor='none', \
                         label = 'GR_M{:.1f}_Rs{:d}'.format \
                                   ( M[m], Rs[rs] ) )
            ax.errorbar( Rs[rs], T_NR   [m,rs], \
                         yerr = T_err_NR[m,rs], \
                         fmt = c[m] + s[2], markerfacecolor='none' , \
                         label = 'NR_M{:.1f}_Rs{:d}'.format \
                                   ( M[m], Rs[rs] ) )

        else:

            ax.plot    ( Rs[rs], T_Mul  [m,rs], c[m] + s[0] )
            ax.errorbar( Rs[rs], T_GR   [m,rs], \
                         yerr = T_err_GR[m,rs], \
                         fmt = c[m] + s[1], markerfacecolor='none' )
            ax.errorbar( Rs[rs], T_NR   [m,rs], \
                         yerr = T_err_NR[m,rs], \
                         fmt = c[m] + s[2], markerfacecolor='none' )

handles, labels = ax.get_legend_handles_labels()
mapping = [ 0, 3, 4, 1, 5, 6, 2, 7, 8 ]
h = []
l = []
for i in range( len( handles ) ):
    h.append( handles[mapping[i]] )
    l.append( labels [mapping[i]] )
ax.legend(h,l,prop={'size':8})

ax.set_xticks( Rs )
ax.set_xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$', labelpad=-0.007 )
ax.set_ylabel( r'$T\ \left[\mathrm{ms}\right]$' )
#plt.show()
plt.savefig( '/home/kkadoogan/fig.PeriodComparison.png', dpi = 300, \
             bbox_inches = 'tight' )
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
