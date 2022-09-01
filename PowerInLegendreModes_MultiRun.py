#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( './Publication.sty' )

if __name__ == '__main__':

    Root = './'

    fig, axs = plt.subplots( 2, 1 )

    suptitle = 'M2.8_Mdot0.3_Rs120'
    IDs = np.array( ['NR2D_M2.8_Mdot0.3_Rs120', \
                     'GR2D_M2.8_Mdot0.3_Rs120' ], str )

    #suptitle = 'GR_M2.8_Mdot0.3'
    #IDs = np.array( ['GR2D_M2.8_Mdot0.3_Rs120', \
    #                 'GR2D_M2.8_Mdot0.3_Rs150', \
    #                 'GR2D_M2.8_Mdot0.3_Rs180' ], str )

    fig.suptitle( suptitle )
    for ID in range( IDs.shape[0] ):

        tS, RsAve, RsMin, RsMax, P0, P1, P2, P3, P4 \
          = np.loadtxt( Root + '.' + IDs[ID] \
                          + '_DivV2_0.80-0.90_PowersInLegendreModes.dat' )

        ind = np.where( tS < 1.5e2 )[0]

        tS    = np.copy( tS   [ind] )
        RsAve = np.copy( RsAve[ind] )
        P1    = np.copy( P1   [ind] )

        if IDs[ID][0:2] == 'GR':
            xscale = np.linspace( 0.8, 1.0, 5 )
            for xs in xscale:
                axs[0].plot( xs*tS, ( RsAve - RsAve[0] ) / RsAve[0], \
                             label = r'$\alpha=${:.2f}'.format( xs ) )
                axs[1].semilogy( xs*tS, P1 )
        else:
            axs[0].plot( tS, ( RsAve - RsAve[0] ) / RsAve[0], \
                         label = 'NR' )
            axs[1].semilogy( tS, P1 )

    axs[0].legend()
    axs[0].xaxis.set_ticklabels([])
#    axs[0].set_ylim( 179.0, 182.0 )
    axs[0].grid()
    axs[1].grid()
    axs[1].set_xlabel( 'Time [ms]' )
    axs[1].set_ylabel( r'$P_{1}$ [cgs]' )
    axs[0].set_ylabel( r'$\left(\left<R_{s}\right>-\left<R_{s,0}\right>\right)/\left<R_{s,0}\right>$' )

    plt.subplots_adjust( hspace = 0.0 )

#    plt.show()
    plt.savefig( 'fig.PowerInLegendreModes_MultiRun_{:}.png'.format \
                 ( suptitle ), dpi = 300 )
    plt.close()

