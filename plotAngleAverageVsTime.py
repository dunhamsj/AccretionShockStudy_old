#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os

from MakeDataFile_new import MakeDataFile, ReadHeader
from UtilitiesModule import GetFileArray

def ComputeAngleAverage( nX, theta, Data, dtheta ):

    iX1 = 310

    uK = 0.5 * np.sum( Data[iX1] * np.sin( theta ) ) * dtheta

    return uK

def PlotAngleAverageVsTime( runID, field ):

    rootDirectory = '/lump/data/accretionShockStudy/'

    plotFileDirectory = rootDirectory + '{:}/'.format( runID )

    UseLogScale = False

    plotFileBaseName = runID + '.plt_'

    ID = '{:}_{:}'.format( runID, field )

    saveFigAs = 'fig.{:}_AAvsTime.png'.format( ID )

    dataFileDirectory = '.{:}/'.format( ID )

    MakeDataFile( field, plotFileDirectory, dataFileDirectory, \
                  plotFileBaseName, 'spherical', Verbose = True )

    plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )

    nSS = plotFileArray.shape[0]
    SSi = 0
    SSf = nSS-1

    uK   = np.empty( (nSS), np.float64 )
    time = np.empty( (nSS), np.float64 )

    for i in range( nSS ):

        iSS = SSi + np.int64( ( SSf - SSi ) / ( nSS - 1 ) * i )

        dataFile = dataFileDirectory + plotFileArray[iSS] + '.dat'

        dataShape, dataUnits, time[i], X1_C, X2_C, X3_C, dX1, dX2, dX3 \
          = ReadHeader( dataFile )

        data = np.loadtxt( dataFile ).reshape( np.int64( dataShape ) )

        nX = X1_C.shape[0]
        uK[iSS] = ComputeAngleAverage( nX, X2_C, data, dX2[0] )

    # Plotting

    # Get mesh data
    dataFile = dataFileDirectory + plotFileArray[0] + '.dat'
    dataShape, dataUnit, t, X1_C, X2_C, X3_C, dX1, dX2, dX3 \
      = ReadHeader( dataFile )

    fig, ax = plt.subplots()

    ax.plot( time, uK, 'k.' )

    ax.set_xlim( time.min(), time.max() )
    ax.set_xlabel( 'Time [ms]' )

    ax.set_ylim( uK.min(), uK.max() )
    ax.set_ylabel( '<{:}> [{:}]'.format( field, dataUnit ) )

    ax.grid()

    #plt.show()
    plt.savefig( saveFigAs, dpi = 300 )

    os.system( 'rm -rf __pycache__ ' )

ID    = 'GR2D_M2.8_Mdot0.3_Rs120'
Field = 'RadialFlux'

PlotAngleAverageVsTime( ID, Field )
