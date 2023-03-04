#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import LogNorm, SymLogNorm

from UtilitiesModule import GetNorm, GetFileArray, GetFileArray_data, \
        ComputeAngleAverage
from MakeDataFile import MakeDataFile, ReadHeader

def MakeMovie2D( ID0, ID1, field0, field1, SSi = -1, SSf = -1, nSS = -1, \
                 movieName = 'mov.TwoFields.mp4' ):

    """
    Generate a movie from data files created in MakeDataFile.py.
    """

    ############# User input #############

    rootDirectory0 = '/lump/data/accretionShockStudy/'
    rootDirectory1 = '/lump/data/accretionShockStudy/'
    plotFileDirectory0 = rootDirectory0 + '{:}/'.format( ID0 )
    plotFileDirectory1 = rootDirectory1 + '{:}/'.format( ID1 )
    dataFileDirectory0 = '.{:}_{:}_MovieData2D/'.format( ID0, field0 )
    dataFileDirectory1 = '.{:}_{:}_MovieData2D/'.format( ID1, field1 )

    cmap0 = 'viridis' # Color scheme for movie
    cmap1 = 'viridis' # Color scheme for movie

    useLogScale0 = True # Use log scale for field?
    useLogScale1 = True # Use log scale for field?

    MovieRunTime = 10.0 # seconds

    zAxisVertical = True # Orient z-axis

    ############# End of user input #############

    plotFileBaseName0 = ID0 + '.plt'
    plotFileBaseName1 = ID1 + '.plt'

    ID0 = '{:}_{:}'.format( ID0, field0 )
    ID1 = '{:}_{:}'.format( ID1, field1 )

    MakeDataFile( field0, plotFileDirectory0, dataFileDirectory0, \
                  plotFileBaseName0, \
                  SSi = SSi, SSf = SSf, nSS = nSS, \
                  verbose = True, \
                  forceChoice = False, OW = True )
    MakeDataFile( field1, plotFileDirectory1, dataFileDirectory1, \
                  plotFileBaseName1, \
                  SSi = SSi, SSf = SSf, nSS = nSS, \
                  verbose = True, \
                  forceChoice = False, OW = True )

    # Get mesh
    #plotFileArray0 = GetFileArray_data( dataFileDirectory0 )
    plotFileArray0 = GetFileArray( plotFileDirectory0, plotFileBaseName0 )
    dataFile0 = dataFileDirectory0 + '{:}'.format( plotFileArray0[0] ) + '.dat'
    dataShape0, dataUnits0, Time0, X1_C0, X2_C0, X3_C0, dX10, dX20, dX30 \
      = ReadHeader( dataFile0 )

    #plotFileArray1 = GetFileArray_data( dataFileDirectory1 )
    plotFileArray1 = GetFileArray( plotFileDirectory1, plotFileBaseName1 )
    dataFile1 = dataFileDirectory1 + '{:}'.format( plotFileArray1[0] ) + '.dat'
    dataShape1, dataUnits1, Time1, X1_C1, X2_C1, X3_C1, dX11, dX21, dX31 \
      = ReadHeader( dataFile1 )

    nR0 = dX10.shape[0]
    nT0 = dX20.shape[0]
    nR1 = dX11.shape[0]
    nT1 = dX21.shape[0]

    if( nSS < 0 ): nSS = min( plotFileArray0.shape[0], plotFileArray1.shape[0] )

    fig      = plt.figure( figsize = (16,9) )
    ax       = fig.add_subplot( 111, polar = True )
    theta0, r0 = np.meshgrid( X2_C0, X1_C0 )
    theta1, r1 = np.meshgrid( 2.0 * np.pi - X2_C1, X1_C1 )

    uK0 = np.empty( (nSS,nR0), np.float64 )
    uK1 = np.empty( (nSS,nR1), np.float64 )

    nodalData0 = np.empty( (nSS,nR0,nT0), np.float64 )
    time0      = np.empty( (nSS)        , np.float64 )
    nodalData1 = np.empty( (nSS,nR1,nT1), np.float64 )
    time1      = np.empty( (nSS)        , np.float64 )

    vmin0 = +np.inf
    vmax0 = -np.inf
    vmin1 = +np.inf
    vmax1 = -np.inf

    for j in range( nSS ):

        dataFile0 = dataFileDirectory0 + plotFileArray0[j] + '.dat'
        dataFile1 = dataFileDirectory1 + plotFileArray1[j] + '.dat'

        dataShape0, dataUnits0, time0[j], \
        X1_C0, X2_C0, X3_C0, dX10, dX20, dX30 \
          = ReadHeader( dataFile0 )
        dataShape1, dataUnits1, time1[j], \
        X1_C1, X2_C1, X3_C1, dX11, dX21, dX31 \
          = ReadHeader( dataFile1 )

        data0 \
          = np.loadtxt( dataFile0, dtype = np.float64 ).reshape( dataShape0 )
        data1 \
          = np.loadtxt( dataFile1, dtype = np.float64 ).reshape( dataShape1 )

        nodalData0[j] = data0
        nodalData1[j] = data1

        nX0 = [ dX10.shape[0], dX20.shape[0], dX30.shape[0] ]
        nX1 = [ dX11.shape[0], dX21.shape[0], dX31.shape[0] ]
        uK0[j] = ComputeAngleAverage( data0, X2_C0, dX20, dX30, nX = nX0 )
        uK1[j] = ComputeAngleAverage( data1, X2_C1, dX21, dX31, nX = nX1 )

        vmin0 = min( vmin0, data0.min() )
        vmax0 = max( vmax0, data0.max() )
        vmin1 = min( vmin1, data1.min() )
        vmax1 = max( vmax1, data1.max() )

    # END for j in range( nSS )

    def f0(t):
        return nodalData0[t]
    def f1(t):
        return nodalData1[t]

    vmin = min( vmin0, vmin1 )
    vmax = max( vmax0, vmax1 )

    #va = abs( min( vmin, -vmax ) )
    #vb = abs( max( vmax, -vmin ) )

    #print( vmin, vmax )
    #vmax = +max( va, vb )
    #vmin = -max( va, vb )
    #print( vmin, vmax )

    Norm0 = GetNorm( useLogScale0, f0(0), vmin = vmin0, vmax = vmax0, \
                    linthresh = 1.0e1 )
    Norm1 = GetNorm( useLogScale1, f1(0), vmin = vmin1, vmax = vmax1, \
                    linthresh = 1.0e1 )

    # Taken from:
    # https://brushingupscience.com/2016/06/21/
    # matplotlib-animations-the-easy-way/
    ax.grid( False )
    cb0axes = fig.add_axes( [ 0.81,  0.1, 0.03, 0.8 ] )
    cb1axes = fig.add_axes( [ 0.185, 0.1, 0.03, 0.8 ] )

    im0 = ax.pcolormesh( theta0, r0, f0(0), \
                         cmap = cmap0, \
                         norm = Norm0, \
                         shading = 'nearest' )
    cbar0 = fig.colorbar( im0, cax = cb0axes )

    im1 = ax.pcolormesh( theta1, r1, f1(0), \
                         cmap = cmap1, \
                         norm = Norm1, \
                         shading = 'nearest' )
    cbar1 = fig.colorbar( im1, cax = cb1axes )

    # Limits on coordinate axes

    rmax = 100.0
    ax.set_thetamin( 0.0 )
    ax.set_thetamax( 360.0 )
    ax.set_rmax( rmax )
    ax.set_theta_direction( -1 )

    if( zAxisVertical ):
        ax.set_theta_zero_location( 'N' ) # z-axis vertical
        time_text0 = ax.text( +0.25 * np.pi / 2.0, rmax * ( 1.0 + 0.1 ), '' )
        time_text1 = ax.text( -0.25 * np.pi / 2.0, rmax * ( 1.0 + 0.1 ), '' )
    else:
        ax.set_theta_zero_location( 'W' ) # z-axis horizontal
        time_text0 = ax.text( 0.9 * np.pi / 2.0, rmax * ( 1.0 + 0.3 ), '' )
        time_text1 = ax.text( 0.3 * np.pi / 2.0, rmax * ( 1.0 + 0.3 ), '' )

    import matplotlib.ticker as ticker
    import matplotlib as mpl

    def UpdateFrame(t):
        im0.set_array( f0(t).flatten() )
        im1.set_array( f1(t).flatten() )
        time_text0.set_text( 'time = {:d} ms'.format( np.int64( time0[t] ) ) )
        time_text1.set_text( 'time = {:d} ms'.format( np.int64( time1[t] ) ) )
        return im0, im1, time_text0, time_text1

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

#field = 'LateralMomentumFluxInRadialDirectionGR'
field = 'PolytropicConstant'
ID1 = 'NR2D_M2.8_Mdot0.3_Rs6.00e1_RPNS2.00e1'
ID2 = 'NR2D_M2.8_Mdot0.5_Rs6.00e1_RPNS2.00e1'

movieName = 'mov.TwoFields.mp4'
MakeMovie2D( ID1, ID2, field, field, movieName = movieName )
