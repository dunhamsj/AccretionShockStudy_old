#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import LogNorm, SymLogNorm

from UtilitiesModule import GetNorm, GetFileArray, ComputeAngleAverage
from MakeDataFile import MakeDataFile, ReadHeader

def MakeMovie2D( ID, field, SSi = -1, SSf = -1, nSS = -1, \
                 movieName = 'OneField.mp4' ):

    """
    Generate a movie from data files created in MakeDataFile.py.
    """

    ############# User input #############

    rootDirectory = '/lump/data/accretionShockStudy/'
    plotFileDirectory = rootDirectory + '{:}/'.format( ID )
    dataFileDirectory = '.{:}_{:}_MovieData2D/'.format( ID, field )

    cmap = 'viridis' # Color scheme for movie

    useLogScale = True # Use log scale for field?

    MovieRunTime = 10.0 # seconds

    zAxisVertical = True # Orient z-axis

    ############# End of user input #############

    plotFileBaseName = ID + '.plt'

    ID = '{:}_{:}'.format( ID, field )

    MakeDataFile( field, plotFileDirectory, dataFileDirectory, \
                  plotFileBaseName, \
                  SSi = SSi, SSf = SSf, nSS = nSS, \
                  verbose = True, \
                  forceChoice = False, OW = True )

    # Get mesh
    plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )[:-1]
    dataFile = dataFileDirectory + '{:}'.format( plotFileArray[0] ) + '.dat'
    dataShape, dataUnits, Time, X1_C, X2_C, X3_C, dX1, dX2, dX3 \
      = ReadHeader( dataFile )

    nR = dX1.shape[0]
    nT = dX2.shape[0]

    fig      = plt.figure( figsize = (16,9) )
    ax       = fig.add_subplot( 111, polar = True )
    theta, r = np.meshgrid( X2_C, X1_C )

    if( nSS < 0 ): nSS = plotFileArray.shape[0]

    uK = np.empty( (nSS,nR), np.float64 )

    nodalData = np.empty( (nSS,nR,nT), np.float64 )
    time      = np.empty( (nSS)      , np.float64 )

    for j in range( nSS ):

        iSS = SSi + np.int64( ( SSf - SSi ) / ( nSS - 1 ) * j )

        dataFile = dataFileDirectory + plotFileArray[iSS] + '.dat'

        dataShape, dataUnits, time[j], X1_C, X2_C, X3_C, dX1, dX2, dX3 \
          = ReadHeader( dataFile )

        data \
          = np.loadtxt( dataFile, dtype = np.float64 ).reshape( dataShape )

        nodalData[j] = data

        uK[j] = ComputeAngleAverage( data, X2_C, dX2, dX3 )

    # END for j in range( nSS )

    def f(t):
        return nodalData[t]

    Norm = GetNorm( useLogScale, f(0) )

    # Taken from:
    # https://brushingupscience.com/2016/06/21/
    # matplotlib-animations-the-easy-way/
    im = ax.pcolormesh( theta, r, f(0), \
                        cmap = cmap, \
                        norm = Norm, \
                        shading = 'nearest' )

    cbar = fig.colorbar( im )

    # Limits on coordinate axes

    rmax = 100.0
    ax.set_thetamin( 0.0 )
    ax.set_thetamax( 180.0 )
    ax.set_rmax( rmax )
    ax.set_theta_direction( -1 )

    if( zAxisVertical ):
        ax.set_theta_zero_location( 'N' ) # z-axis vertical
        time_text = plt.text( 0.125 * np.pi / 2.0, rmax * ( 1.0 + 0.1 ), '' )
    else:
        ax.set_theta_zero_location( 'W' ) # z-axis horizontal
        time_text = plt.text( 0.9 * np.pi / 2.0, rmax * ( 1.0 + 0.3 ), '' )

    import matplotlib.ticker as ticker
    import matplotlib as mpl

    def UpdateFrame(t):
        im.set_array( f(t).flatten() )
        time_text.set_text( 'time = {:d} ms'.format( np.int64( time[t] ) ) )
        return im, time_text

    # Call the animator

    print( 'Making movie...' )

    fps = max( np.int64( nSS / MovieRunTime ), 1 )

    anim \
      = animation.FuncAnimation \
          ( fig, UpdateFrame, frames = nSS, blit = True )

    anim.save( movieName, fps = fps )
    plt.close()

    import os
    os.system( 'rm -rf __pycache__ ' )

field = 'NonRadialKineticEnergyDensityGR'
ID = 'GR2D_M2.8_Mdot0.3_Rs6.00e1_RPNS2.00e1'

movieName \
  = 'mov.{:}_{:}.mp4'.format( ID, field )
MakeMovie2D( ID, field, movieName = movieName )
