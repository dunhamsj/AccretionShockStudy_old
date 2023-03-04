#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

rel    = [ 'NR', 'GR' ]
M      = [ '2.8' ]
Mdot   = [ '0.3' ]
Rs     = [ '6.00e1', '7.50e1', '9.00e1' ]
suffix = '_RPNS2.00e1'

fig, axs = plt.subplots( 2, 3, figsize = (10,6) )

for r in range( len( rel ) ):
    for m in M:
        for mdot in Mdot:
            for rs in range( len( Rs ) ):

                ID = '{:}2D_M{:}_Mdot{:}_Rs{:}{:}' \
                     .format( rel[r], m, mdot, Rs[rs], suffix )

                dataFileName = '{:}_ShockRadiusVsTime.dat'.format( ID )

                Time, RsAve, RsMin, RsMax = np.loadtxt( dataFileName )

                dr = ( RsAve - RsAve[0] ) / RsAve[0]

                if r == 0:
                    axs[r,rs].set_title( r'$R_{{s}}={:}\ \mathrm{{km}}$' \
                                         .format( Rs[rs] ) )

                if rel[r] == 'NR':

                    axs[r,rs].plot( Time, dr, \
                                    'k-', label = 'RsAve (NR)' )
                else:

                    axs[r,rs].plot( Time, dr, \
                                    'k-', label = 'RsAve (GR)' )

                axs[r,rs].grid()

for i in range( axs.shape[0] ):
    axs[i,0].set_ylabel \
               ( r'$\left(R_{s}\left(t\right)-R_{s}\left(0\right)\right)$' \
               + r'$/R_{s}\left(0\right)$', labelpad = +0.2 )
for j in range( axs.shape[1] ):
    axs[1,j].set_xlabel( 'Time [ms]' )

fig.suptitle( r'$M=2.8\ \mathrm{M}_{\odot}$' )
axs[0,0].legend()
axs[1,0].legend()
plt.subplots_adjust( hspace = 0.0 )
plt.show()
#plt.savefig( 'fig.ShockRadiusVsTime_LateStage.png', dpi = 300 )
plt.close()
