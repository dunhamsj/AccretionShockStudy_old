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

c = [ 'r', 'b', 'm' ]
s = [ 'x', '.', 's' ]

NX = 5.0
NY = NX * 3/4
fig, ax = plt.subplots( 1, 1, figsize = (NX,NY) )

for rs in range( Rs.shape[0] ):
    for m in range( M.shape[0] ):

        if rs == 0:

            ax.errorbar( Rs[rs], G_GR   [m,rs], \
                         yerr = G_err_GR[m,rs], \
                         fmt = c[m] + s[1], markerfacecolor='none', \
                         label = 'GR_M{:.1f}_Rs{:d}'.format \
                                   ( M[m], Rs[rs] ) )
            ax.errorbar( Rs[rs], G_NR   [m,rs], \
                         yerr = G_err_NR[m,rs], \
                         fmt = c[m] + s[2], markerfacecolor='none' , \
                         label = 'NR_M{:.1f}_Rs{:d}'.format \
                                   ( M[m], Rs[rs] ) )

        else:

            ax.errorbar( Rs[rs], G_GR   [m,rs], \
                         yerr = G_err_GR[m,rs], \
                         fmt = c[m] + s[1], markerfacecolor='none' )
            ax.errorbar( Rs[rs], G_NR   [m,rs], \
                         yerr = G_err_NR[m,rs], \
                         fmt = c[m] + s[2], markerfacecolor='none' )

ax.legend()

ax.set_xticks( Rs )
ax.set_xlabel( r'$R_{s}\ \left[\mathrm{km}\right]$', labelpad=-0.007 )
ax.set_ylabel( r'$\omega_{r}\ \left[\mathrm{1/ms}\right]$' )
#plt.show()
plt.savefig( '/home/kkadoogan/fig.GrowthRateComparison.png', dpi = 300, \
             bbox_inches = 'tight' )
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
