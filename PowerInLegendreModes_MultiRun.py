#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( './Publication.sty' )

if __name__ == '__main__':

    Root = './'

    IDs = np.array( ['NR2D_M0.14_Mdot0.03_Rs180_PA1.00e-06_nX640x064', \
                     'GR2D_M0.14_Mdot0.03_Rs180_PA1.00e-06_nX640x064' ], str )

    fig, axs = plt.subplots( 2, 1 )

    fig.suptitle( IDs[0][5:] )

    for ID in range( IDs.shape[0] ):

        tS, RsAve, RsMin, RsMax, P0, P1, P2, P3, P4 \
          = np.loadtxt( Root + '.' + IDs[ID] \
                          + '_DivV2_0.80-0.90_PowersInLegendreModes.dat' )

        ind = np.where( tS < 3.0e2 )[0]

        tS    = np.copy( tS   [ind] )
        RsAve = np.copy( RsAve[ind] )
        P1    = np.copy( P1   [ind] )

        axs[0].plot( tS, RsAve, label = IDs[ID][0:2] )
        axs[1].semilogy( tS, P1 )

    axs[0].legend()
    axs[0].xaxis.set_ticklabels([])
    axs[0].set_ylim( 179.0, 182.0 )
    axs[0].grid()
    axs[1].grid()
    axs[1].set_xlabel( 'Time [ms]' )
    axs[1].set_ylabel( r'$P_{1}$ [cgs]' )
    axs[0].set_ylabel( r'$\left<R_{s}\right>$ [km]' )

    plt.subplots_adjust( hspace = 0.0 )

#    plt.show()
    plt.savefig( 'fig.PowerInLegendreModes_MultiRun.png', dpi = 300 )
    plt.close()

