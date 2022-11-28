#!/usr/bin/env python3

import yt
import numpy as np
import matplotlib.pyplot as plt
import os

from UtilitiesModule import Overwrite, GetData, GetFileArray

yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

def MakeLineOutPlot( plotFileDirectory, plotFileBaseName, entropyThreshold ):

    plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )

    plotFile = plotFileDirectory + plotFileArray[-2]
    time, data, dmy1, r, dmy3, dmy4, dmy5, dmy6, dmy7, dmy8 \
      = GetData( plotFile, 'PolytropicConstant' )

    fig, ax = plt.subplots()
    ax.semilogy( r, data[:,0,0], 'k-' )
    ax.text( 0.5, 0.7, 'Time = {:.3e}'.format( time ), \
             transform = ax.transAxes )
    ax.axhline( entropyThreshold, label = 'Shock radius cut-off' )

    plt.legend()
#    plt.savefig( 'entropyThresholdCheck.png', dpi = 300 )
    plt.show()
    plt.close()

    return


def MakeDataFile \
      ( plotFileDirectory, plotFileBaseName, dataFileName, \
        entropyThreshold, markEvery = 1, forceChoice = False, OW = True ):

    OW = Overwrite( dataFileName, ForceChoice = forceChoice, OW = OW )

    if not OW: return

    print( '\n  Creating {:}'.format( dataFileName ) )
    print(   '  --------' )

    plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )

    plotFileArray = plotFileArray[0::markEvery]

    # Just to get number of elements
    plotFile = plotFileDirectory + plotFileArray[0]
    dmy0, dmy1, dmy2, dmy3, dmy4, dmy5, dmy6, dmy7, dmy8, nX \
      = GetData( plotFile, 'PolytropicConstant' )

    nSS = plotFileArray.shape[0]

    Data = np.empty( (nSS,nX[0],nX[1],nX[2]), np.float64 )

    Time = np.empty( nSS, np.float64 )

    Volume = np.zeros( nSS, np.float64 )
    RsMin  = np.empty( nSS, np.float64 )
    RsMax  = np.empty( nSS, np.float64 )

    for i in range( nSS ):

        print( '    {:}/{:}'.format( i,nSS) )

        plotFile = plotFileDirectory + plotFileArray[i]

        Time[i], Data[i], dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFile, 'PolytropicConstant' )

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

    rootDirectory = '/lump/data/accretionShockStudy/'
    #rootDirectory = '/home/dunhamsj/Work/thornado_GW/SandBox/AMReX/Applications/StandingAccretionShock_Relativistic/'
    ID = 'GR2D_M2.8_Mdot0.3_Rs180'
    plotFileDirectory = rootDirectory + ID + '/'
    plotFileBaseName = ID + '.plt'
    entropyThreshold = 1.0e15

#    MakeLineOutPlot( plotFileDirectory, plotFileBaseName, entropyThreshold )

    dataFileName = '{:}_ShockRadiusVsTime.dat'.format( ID )
    forceChoice = False
    OW = False
    MakeDataFile \
      ( plotFileDirectory, plotFileBaseName, dataFileName, \
        entropyThreshold, markEvery = 1, forceChoice = forceChoice, OW = OW )

    Time, RsAve, RsMin, RsMax = np.loadtxt( dataFileName )

    plt.plot( Time, ( RsAve - RsAve[0] ) / RsAve[0], 'k-', label = 'RsAve' )
#    plt.plot( Time, RsMin, 'r-', label = 'RsMin' )
#    plt.plot( Time, RsMax, 'b-', label = 'RsMax' )
    #plt.fill_between( Time, RsMin, RsMax )
    plt.legend()
#    plt.savefig( 'fig.{:}_ShockRadiusVsTime.png'.format( ID ), dpi = 300 )
    plt.show()
    plt.close()
    import os
    os.system( 'rm -rf __pycache__ ' )
