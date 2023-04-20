#!/usr/bin/env python3

import yt
import numpy as np
import matplotlib.pyplot as plt
import os
plt.style.use( 'publication.sty' )

from UtilitiesModule import Overwrite, GetData, GetFileArray

yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

def MakeLineOutPlot( plotfileDirectory, plotfileBaseName, entropyThreshold ):

    plotfileArray = GetFileArray( plotfileDirectory, plotfileBaseName )

    data, DataUnit, r, X2_C, X3_C, dX1, dX2, dX3, xL, xH, nX, time \
      = GetData( plotfileDirectory, plotfileBaseName, 'PolytropicConstant', \
                 'spherical', True, argv = ['a','0'], \
                 ReturnTime = True, ReturnMesh = True, Verbose = True )

    fig2, ax2 = plt.subplots()
    ax2.semilogy( r[:,0,0], data[:,0,0], 'k-' )
    ax2.text( 0.5, 0.7, 'Time = {:.3e}'.format( time ), \
              transform = ax2.transAxes )
    ax2.axhline( entropyThreshold, label = 'Shock radius cut-off' )

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

    print( '\n plotfileDirectory: ', plotfileDirectory )
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
    rootDirectory = '/lump/data/accretionShockStudy/newData/2D/'

#    M    = [ '1.4' ]
#    Rs   = [ '1.20e2', '1.50e2', '1.75e2' ]
#    xL   = [ '040' ]

    M    = [ '2.8' ]
    Rs   = [ '6.00e1', '7.00e1', '8.75e1' ]
    xL   = [ '020' ]

    rel  = [ 'NR', 'GR', ]
    Mdot = [ '0.3' ]
    suffix = [ '' ]
    xH   = [ 1.5 * x for x in np.float64( Rs ) ]

    fig, ax = plt.subplots( 1, 1 )

    # colorblind-friendly palette: https://gist.github.com/thriveth/8560036
    color = ['#377eb8', '#ff7f00', '#4daf4a', \
             '#f781bf', '#a65628', '#984ea3', \
             '#999999', '#e41a1c', '#dede00']

    for i in range( len( rel ) ):
        for j in range( len( Rs ) ):

            ID = '{:}2D_M{:}_Rpns{:}_Rs{:}{:}' \
                 .format( rel[i], M[0], xL[0], Rs[j], suffix[0] )

            plotfileDirectory = rootDirectory + ID + '/'
            plotfileBaseName = ID + '.plt'
            entropyThreshold = 1.0e15

            #MakeLineOutPlot \
            #  ( plotfileDirectory, plotfileBaseName, entropyThreshold )

            dataFileName = '.{:}_ShockRadiusVsTime.dat'.format( ID )
            forceChoice = True
            OW = False
            MakeDataFile \
              ( plotfileDirectory, plotfileBaseName, dataFileName, \
                entropyThreshold, markEvery = 1, forceChoice = forceChoice, \
                OW = OW )

            Time, RsAve, RsMin, RsMax = np.loadtxt( dataFileName )
            Time  = np.copy( Time [:-1] )
            RsAve = np.copy( RsAve[:-1] )
            RsMin = np.copy( RsMin[:-1] )
            RsMax = np.copy( RsMax[:-1] )

            ind = -1
            ind = np.where( RsMax > 1.45 * np.float64( Rs[j] ) )[0]
            if not len( ind ) == 0:
                ind = np.copy( ind[0] )
            else:
                ind = -1
            print( ID, ind, Time[ind] )

            c = i
            ax.plot( Time[0:ind], RsAve[0:ind] / RsAve[0], \
                    c = color[c], ls = '-' , \
                    label = '{:}'.format( rel[i] ) )
            ax.plot( Time[0:ind], RsMin[0:ind] / RsAve[0], \
                     c = color[c], ls = '--', label = 'min' )
            ax.plot( Time[0:ind], RsMax[0:ind] / RsAve[0], \
                     c = color[c], ls = ':' , label = 'max' )

    ax.set_xlabel( 'Time [ms]' )
    ax.set_ylabel( r'$R_{\mathrm{S}}/R_{\mathrm{S}}\left(0\right)$', \
                   labelpad = +10.0 )
    #ax.axhline( 1.01 )
    #ax.axhline( 0.99 )
    ax.grid()
    ax.legend()

    #plt.savefig( 'fig.{:}_ShockRadiusVsTime.png'.format( ID ), dpi = 300 )
    plt.show()
    plt.close()

    import os
    os.system( 'rm -rf __pycache__ ' )
