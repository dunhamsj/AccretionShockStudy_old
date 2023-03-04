#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

def PlotRelaxationVsTime \
      ( ax, Time, Data, field, ID, label = '' ):

    return

if __name__ == '__main__':

    nX = [ 128, 256, 512, 1024 ]

    UseLogScale = True

    ID = 'GR1D_M2.8_Mdot0.3_Rs7.50e1_RPNS2.00e1'

    SaveFileAs = 'fig.Relaxation_{:}.png'.format( ID )

#    Root = '/home/dunhamsj/AccretionShockData/'
    Root = '/lump/data/accretionShockStudy/'

    fig, axs = plt.subplots( 3, 1 )

    fig.suptitle( r'$\texttt{{{:}}}$'.format( ID ) )

    for n in nX:

        dataFileName = '.{:}_Relaxation_{:}_nX{:}.dat' \
                       .format( ID, 'PF_D', str( n ).zfill(4) )
        Time, D = np.loadtxt( dataFileName )

        dataFileName = '.{:}_Relaxation_{:}_nX{:}.dat' \
                       .format( ID, 'PF_V1', str( n ).zfill(4) )
        Time, V = np.loadtxt( dataFileName )

        dataFileName = '.{:}_Relaxation_{:}_nX{:}.dat' \
                       .format( ID, 'AF_P', str( n ).zfill(4) )
        Time, P = np.loadtxt( dataFileName )

        axs[0].plot( Time, np.abs( D ), '.', \
                 markersize = 2.0, markevery = 1, label = n )
        axs[1].plot( Time, np.abs( V ), '.', \
                 markersize = 2.0, markevery = 1 )
        axs[2].plot( Time, np.abs( P ), '.', \
                 markersize = 2.0, markevery = 1 )

    for i in range( axs.shape[0] ):
        axs[i].grid()
        axs[i].set_yscale( 'log' )

    axs[0].legend()
    axs[-1].set_xlabel( 'Time [ms]' )
    axs[0].set_ylabel( r'$\mathrm{max}\left(\left|\dot{\rho}/\rho\right|\right)$' )
    axs[1].set_ylabel( r'$\mathrm{max}\left(\left|\dot{v}/v\right|\right)$' )
    axs[2].set_ylabel( r'$\mathrm{max}\left(\left|\dot{p}/p\right|\right)$' )

#    plt.show()
    plt.savefig( 'fig.Relaxation_{:}.png'.format( ID ), dpi = 300 )

    plt.close()

    import os
    os.system( 'rm -rf __pycache__ ' )
