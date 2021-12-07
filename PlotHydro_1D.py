#!/usr/bin/env python3

import matplotlib.pyplot as plt
from sys import argv

from UtilitiesModule import GetData, GetNorm

def PlotHydro_1D( ID ):

    fig, axs = plt.subplots( 4, 2, figsize = (12,8) )

    RootDirectory \
      = '/Users/dunhamsj/Research/Data/AccretionShockParameterStudy/'

    for i in range( len( ID ) ):

        DataDirectory = RootDirectory + ID[i] + '/'
        PlotFileBaseName = ID[i] + '.plt'

        D, DataUnit, r, theta, Time, xL, xU \
          = GetData( DataDirectory, PlotFileBaseName, argv, 'PF_D' )

        V, DataUnit, r, theta, Time, xL, xU \
          = GetData( DataDirectory, PlotFileBaseName, argv, 'PF_V1' )

        P, DataUnit, r, theta, Time, xL, xU \
          = GetData( DataDirectory, PlotFileBaseName, argv, 'AF_P' )

        C, DataUnit, r, theta, Time, xL, xU \
          = GetData( DataDirectory, PlotFileBaseName, argv, 'AF_Cs' )

        iX2 = [ 0 ]

        if ID[i][0] == 'G':

            p, DataUnit, r, theta, Time, xL, xU \
              = GetData( DataDirectory, PlotFileBaseName, argv, 'GF_Psi' )

            a, DataUnit, r, theta, Time, xL, xU \
              = GetData( DataDirectory, PlotFileBaseName, argv, 'GF_Alpha' )

            h, DataUnit, r, theta, Time, xL, xU \
              = GetData( DataDirectory, PlotFileBaseName, argv, \
                         'SpecificEnthalpy' )

            W, DataUnit, r, theta, Time, xL, xU \
              = GetData( DataDirectory, PlotFileBaseName, argv, \
                         'LorentzFactor' )

            axs[0,1].plot( r[:,iX2], p[:,iX2] )
            axs[1,1].plot( r[:,iX2], a[:,iX2] )
            axs[2,1].plot( r[:,iX2], h[:,iX2] )
            axs[3,1].plot( r[:,iX2], W[:,iX2] )

        for i in range( len( iX2 ) ):

            if i == 0:

                axs[0,0].semilogy( r[:,iX2], D[:,iX2], \
                                   label = ID[i] )

            else:

                axs[0,0].semilogy( r[:,iX2], D[:,iX2] )

        axs[1,0].semilogy( r[:,iX2], P[:,iX2] )
        axs[2,0].plot    ( r[:,iX2], V[:,iX2] )
        axs[3,0].plot    ( r[:,iX2], C[:,iX2] )

    axs[0,0].set_ylabel( r'$\rho\,\left[\mathrm{g\,cm}^{-3}\right]$' )
    axs[1,0].set_ylabel( r'$p\,\left[\mathrm{erg\,cm}^{-3}\right]$' )
    axs[2,0].set_ylabel( r'$v\,\left[\mathrm{km\,s}^{-1}\right]$' )
    axs[3,0].set_ylabel( r'$c_{s}\,\left[\mathrm{km\,s}^{-1}\right]$' )

    axs[0,1].set_ylabel( r'$\psi$' )
    axs[1,1].set_ylabel( r'$\alpha$' )
    axs[2,1].set_ylabel( r'$h/c^{2}$' )
    axs[3,1].set_ylabel( r'$W$' )

    axs[0,1].yaxis.tick_right()
    axs[0,1].yaxis.set_label_position("right")
    axs[1,1].yaxis.tick_right()
    axs[1,1].yaxis.set_label_position("right")
    axs[2,1].yaxis.tick_right()
    axs[2,1].yaxis.set_label_position("right")
    axs[3,1].yaxis.tick_right()
    axs[3,1].yaxis.set_label_position("right")

    axs[0,0].xaxis.set_visible( False )
    axs[1,0].xaxis.set_visible( False )
    axs[2,0].xaxis.set_visible( False )
    axs[3,0].set_xlabel( r'Radial Coordinate [km]' )

    axs[0,1].xaxis.set_visible( False )
    axs[1,1].xaxis.set_visible( False )
    axs[2,1].xaxis.set_visible( False )
    axs[3,1].set_xlabel( r'Radial Coordinate [km]' )

    axs[0,0].legend()

    plt.subplots_adjust( hspace = 0.0, wspace = 0.0 )

    #plt.savefig( 'Hydro_1D.png', dpi = 300 )

    plt.show()
    plt.close()

ID = [ 'GR2D_M1.4_Mdot0.3_Rs180_PA0.040_nX320x064', \
       'GR2D_M1.4_Mdot0.3_Rs180_PA0.040_nX320x128' ]

PlotHydro_1D( ID )

import os
os.system( 'rm -rf __pycache__ ' )
