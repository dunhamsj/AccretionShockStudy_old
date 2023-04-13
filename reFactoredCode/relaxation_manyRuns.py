#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

def PlotRelaxationVsTime \
      ( ax, Time, Data, field, ID, label = '' ):

    return

if __name__ == '__main__':

    UseLogScale = True

    ID = [ 'NR1D_M1.4_Rpns040_Rs1.80e2_ST1.00e-10', \
           'NR1D_M1.4_Rpns040_Rs180_Mdot0.3' ]

    ST = [ 'ST = 1.00e-10', 'ST = 1.00e-06' ]
    SaveFileAs = 'fig.Relaxation_{:}.png'.format( ID )

    fig, ax = plt.subplots( 1, 1 )

    for i in range( len( ID ) ):

        dataFileName = '.{:}_Relaxation_{:}_nX0460.dat' \
                       .format( ID[i], 'PF_D' )
        Time, D = np.loadtxt( dataFileName )

        ax.plot( Time, np.abs( D ), '.', \
                 markersize = 3.0-2*i, markevery = 1, label = ST[i] )

        ax.grid()
        ax.set_yscale( 'log' )

    ax.legend()
    ax.set_xlabel( 'Time [ms]' )
    ax.set_ylabel( r'$\mathrm{max}\left(\left|\dot{\rho}/\rho\right|\right)$' )

#    plt.show()
    plt.savefig( 'fig.Relaxation_{:}.png'.format( 'ST' ), dpi = 300 )

    plt.close()

    import os
    os.system( 'rm -rf __pycache__ ' )
