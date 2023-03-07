#!/usr/bin/env python3

import yt
import numpy as np
import matplotlib.pyplot as plt
import os

from UtilitiesModule import Overwrite, GetData, GetFileArray

yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

def MakeLineOutPlot( plotfileDirectory, plotfileBaseName, entropyThreshold ):

    plotfileArray = GetFileArray( plotfileDirectory, plotfileBaseName )

    data, DataUnit, r, X2_C, X3_C, dX1, dX2, dX3, xL, xH, nX, time \
      = GetData( plotfileDirectory, plotfileBaseName, 'PolytropicConstant', \
                 'spherical', True, \
                 ReturnTime = True, ReturnMesh = True, Verbose = True )

    fig, ax = plt.subplots()
    ax.semilogy( r[:,0,0], data[:,0,0], 'k-' )
    ax.text( 0.5, 0.7, 'Time = {:.3e}'.format( time ), \
             transform = ax.transAxes )
    ax.axhline( entropyThreshold, label = 'Shock radius cut-off' )

    plt.legend()
#    plt.savefig( 'entropyThresholdCheck.png', dpi = 300 )
    plt.show()
    plt.close()

    return


def MakeDataFile \
      ( plotfileDirectory, plotfileBaseName, dataFileName, \
        entropyThreshold, markEvery = 1, forceChoice = False, OW = True ):

    OW = Overwrite( dataFileName, ForceChoice = forceChoice, OW = OW )

    if not OW: return

    print( '\n  Creating {:}'.format( dataFileName ) )
    print(   '  --------' )

    plotfileArray = GetFileArray( plotfileDirectory, plotfileBaseName )

    plotfileArray = plotfileArray[0::markEvery]

    data, DataUnit, r, X2_C, X3_C, dX1, dX2, dX3, xL, xH, nX, time \
      = GetData( plotfileDirectory, plotfileBaseName, 'PolytropicConstant', \
                 'spherical', True, \
                 ReturnTime = True, ReturnMesh = True, Verbose = False )

    nSS = plotfileArray.shape[0]

    Data = np.empty( (nSS,nX[0],nX[1],nX[2]), np.float64 )

    Time = np.empty( nSS, np.float64 )

    Volume = np.zeros( nSS, np.float64 )
    RsMin  = np.empty( nSS, np.float64 )
    RsMax  = np.empty( nSS, np.float64 )

    for i in range( nSS ):

        print( '    {:}/{:}'.format( i,nSS) )

        plotfile = plotfileDirectory + plotfileArray[i]

        Data[i], dataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xH, nX, Time[i] \
          = GetData( plotfileDirectory, plotfileBaseName, \
                     'PolytropicConstant', \
                     'spherical', True, argv = ['a',plotfile[-8:]], \
                     ReturnTime = True, ReturnMesh = True, Verbose = False )

        X1 = np.copy( X1[:,0,0] )
        X2 = np.copy( X2[0,:,0] )
        dX1 = np.copy( dX1[:,0,0] )
        dX2 = np.copy( dX2[0,:,0] )

        ind = np.where( Data[i] < entropyThreshold )[0]
        RsMin[i] = X1[ind.min()]

        ind = np.where( Data[i] > entropyThreshold )[0]
        RsMax[i] = X1[ind.max()]

        for iX1 in range( nX[0] ):
            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):

                    if( Data[i,iX1,iX2,iX3] > entropyThreshold ):

                        if( nX[1] > 1 ):

                            Volume[i] \
                              += 2.0 * np.pi \
                                   * 1.0 / 3.0 * ( ( X1[iX1] + dX1[iX1] )**3 \
                                                       - X1[iX1]**3 ) \
                                   * ( np.cos( X2[iX2] ) \
                                         - np.cos( X2[iX2] + dX2[iX2] ) )

                        else:

                            Volume[i] \
                              += 4.0 * np.pi \
                                   * 1.0 / 3.0 * ( ( X1[iX1] + dX1[iX1] )**3 \
                                                       - X1[iX1]**3 )

        Rs = ( Volume[i] / ( 4.0 / 3.0 * np.pi ) )**( 1.0 / 3.0 )
        if RsMin[i] >= Rs: RsMin[i] = Rs
        if RsMax[i] <= Rs: RsMax[i] = Rs

    # 4/3 * pi * Rs^3 = Volume
    # ==> Rs = ( Volume / ( 4/3 * pi ) )^( 1 / 3 )
    RsAve = ( Volume / ( 4.0 / 3.0 * np.pi ) )**( 1.0 / 3.0 )

    np.savetxt( dataFileName, \
                np.vstack( ( Time, RsAve, RsMin, RsMax ) ) )

    return


if __name__ == "__main__":

    #rootDirectory = '/lump/data/accretionShockStudy/'
    rootDirectory = '/lump/data/accretionShockStudy/newRuns/'

    rel  = [ 'NR' ]
    M    = [ '1.4' ]
    Mdot = [ '0.3' ]
    Rs   = [ '150' ]
    nX   = [ '640' ]

    fig, ax = plt.subplots( 1, 1 )

    ID = 'NR2D_M1.4_Rpns040_Rs150_Mdot0.3'
    ax.set_title( ID )

    for nx in nX:

#        ID = IDD + '.nX{:}'.format( nx )
        plotfileDirectory = rootDirectory + ID + '/'
        plotfileBaseName = ID + '.plt'
        entropyThreshold = 1.0e15

#        MakeLineOutPlot \
#          ( plotfileDirectory, plotfileBaseName, entropyThreshold )

        dataFileName = '{:}_ShockRadiusVsTime.dat'.format( ID )
        forceChoice = False
        OW = False
        MakeDataFile \
          ( plotfileDirectory, plotfileBaseName, dataFileName, \
            entropyThreshold, markEvery = 10, forceChoice = forceChoice, \
            OW = OW )

        Time, RsAve, RsMin, RsMax = np.loadtxt( dataFileName )

        dr = ( 3.60e2 - 4.00e1 ) / np.float64( nx )

        lab = 'dr = {:.2f} km'.format( dr )
        ax.plot( Time, ( RsAve - RsAve[0] ) / RsAve[0], label = lab )

    ax.set_xlabel( 'Time [ms]' )
    ax.set_ylabel( r'$\left(R_{s}\left(t\right)-R_{s}\left(0\right)\right)/R_{s}\left(0\right)$', labelpad = -0.1 )
    ax.grid()
    ax.legend()

    #plt.savefig( 'fig.{:}_ShockRadiusVsTime.png'.format( ID ), dpi = 300 )
    plt.show()
    plt.close()

    import os
    os.system( 'rm -rf __pycache__ ' )
