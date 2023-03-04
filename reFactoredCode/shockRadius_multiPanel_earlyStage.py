#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

rel    = [ 'NR', 'GR' ]
M      = [ '1.4', '2.0', '2.8' ]
Mdot   = [ '0.3' ]
Rs     = [ '120', '150', '180' ]
suffix = ''

fig, axs = plt.subplots( 3, 3, figsize = (10,6) )

for m in range( len( M ) ):
    for mdot in range( len( Mdot ) ):
        for rs in range( len( Rs ) ):

            ID_NR = '{:}2D_M{:}_Mdot{:}_Rs{:}{:}' \
                    .format( rel[0], M[m], Mdot[mdot], Rs[rs], suffix )
            ID_GR = '{:}2D_M{:}_Mdot{:}_Rs{:}{:}' \
                    .format( rel[1], M[m], Mdot[mdot], Rs[rs], suffix )

            dataFileName_NR = '{:}_ShockRadiusVsTime.dat'.format( ID_NR )
            dataFileName_GR = '{:}_ShockRadiusVsTime.dat'.format( ID_GR )

            Time_NR, RsAve_NR, RsMin_NR, RsMax_NR \
              = np.loadtxt( dataFileName_NR )
            Time_GR, RsAve_GR, RsMin_GR, RsMax_GR \
              = np.loadtxt( dataFileName_GR )

            ind = min( np.argmax( Time_NR ), np.argmax( Time_GR ) )

            dr_NR = ( RsAve_NR - RsAve_NR[0] ) / RsAve_NR[0]
            dr_GR = ( RsAve_GR - RsAve_GR[0] ) / RsAve_GR[0]

            axs[m,rs].plot( Time_NR[:ind], dr_NR[:ind], label = 'NR' )
            axs[m,rs].plot( Time_GR[:ind], dr_GR[:ind], label = 'GR' )

            label = r'$\texttt{{M{:}_Rs{:}}}$'.format( M[m], Rs[rs] )
            axs[m,rs].text( 0.1, 0.2, label, transform = axs[m,rs].transAxes, \
                            fontsize = 12 )

            axs[m,rs].grid()
            axs[m,rs].set_xlim( 0.0, 150.0 )
            if( m < axs.shape[0]-1 ):
                axs[m,rs].set_xticklabels( '' )

fig.supylabel( r'$\left(R_{s}\left(t\right)-R_{s}\left(0\right)\right)$' \
               + r'$/R_{s}\left(0\right)$', x = 0.05 )
fig.supxlabel( 'Time [ms]', y = 0.04 )

axs[0,0].legend()
plt.subplots_adjust( hspace = 0.0 )
plt.show()
#plt.savefig( 'fig.ShockRadiusVsTime_EarlyStage.png', dpi = 300 )
plt.close()
