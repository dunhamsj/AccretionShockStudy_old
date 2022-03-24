#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import LogNorm, SymLogNorm
import os

from MakeDataFile import MakeDataFile
from UtilitiesModule import GetNorm

def ComputeAngleAverage( nX, theta, Data, dt ):

    uK = np.empty( (nX[0]), np.float64 )

    for iX1 in range( uK.shape[0] ):

        uK[iX1] = 0.5 * np.sum( Data[iX1] * np.sin( theta ) ) * dt

    return uK

def MakeMovie2D( ID1, ID2, Field1, Field2, MovieName = 'TwoFields.mp4' ):

    """
    Generate a movie from data files created in MakeDataFile.py.
    """

    ############# User input #############

    Root = '/home/dunhamsj/AccretionShockData'
    DataDirectory1 = Root + '/{:}/'.format( ID1 )
    DataDirectory2 = Root + '/{:}/'.format( ID2 )

    cmap1 = 'viridis' # Color scheme for movie
    cmap2 = 'RdBu' # Color scheme for movie

    UseLogScale1 = True # Use log scale for field?
    UseLogScale2 = True # Use log scale for field?

    MovieRunTime = 10.0 # seconds

    UseCustomTicks1 = False # Define limits for colorbar
    UseCustomTicks2 = False # Define limits for colorbar

    zAxisVertical = True # Orient z-axis

    ############# End of user input #############

    PlotFileBaseName1 = ID1 + '.plt'
    PlotFileBaseName2 = ID2 + '.plt'

    ID1 = '{:}_{:}'.format( ID1, Field1 )
    ID2 = '{:}_{:}'.format( ID2, Field2 )

    DataFileName1 = '.{:}.dat'.format( ID1 )
    DataFileName2 = '.{:}.dat'.format( ID2 )

    xL1, xU1, nX1, FileArray1 \
      = MakeDataFile( Field1, DataDirectory1, DataFileName1, \
                      PlotFileBaseName1)#, SSi = 0, SSf = 600 )

    xL2, xU2, nX2, FileArray2 \
      = MakeDataFile( Field2, DataDirectory2, DataFileName2, \
                      PlotFileBaseName2)#, SSi = 0, SSf = 600 )

    xL1 = xL1
    xU1 = xU1
    dr1 = ( xU1[0] - xL1[0] ) / np.float64( nX1[0] )
    dt1 = ( xU1[1] - xL1[1] ) / np.float64( nX1[1] )

    xL2 = xL2
    xU2 = xU2
    dr2 = ( xU2[0] - xL2[0] ) / np.float64( nX2[0] )
    dt2 = ( xU2[1] - xL2[1] ) / np.float64( nX2[1] )

    f1 = open( DataFileName1 )
    header1 = f1.readline()[16:-2]
    Field1Unit = f1.readline()[9:-1]
    DataShape1 = tuple( [ np.int64( dim ) for dim in header1.split( ',' ) ] )

    f2 = open( DataFileName2 )
    header2 = f2.readline()[16:-2]
    Field2Unit = f2.readline()[9:-1]
    DataShape2 = tuple( [ np.int64( dim ) for dim in header2.split( ',' ) ] )

    Time1 = list( [ np.float64( t ) for t in f1.readline()[12:-2].split(' ') ] )
    Time1 = np.array( Time1 )

    Time2 = list( [ np.float64( t ) for t in f2.readline()[12:-2].split(' ') ] )
    Time2 = np.array( Time2 )

    Data1 \
      = np.loadtxt( DataFileName1, dtype = np.float64 ).reshape( DataShape1 )
    Data2 \
      = np.loadtxt( DataFileName2, dtype = np.float64 ).reshape( DataShape2 )

    fig        = plt.figure( figsize = (16,9) )
    ax         = fig.add_subplot( 111, polar = True )
    X11        = np.linspace( xL1[0] + dr1 / 2.0, xU1[0] - dr1 / 2.0, nX1[0] )
    X21        = np.linspace( xL1[1] + dt1 / 2.0, xU1[1] - dt1 / 2.0, nX1[1] )
    X12        = np.linspace( xL2[0] + dr2 / 2.0, xU2[0] - dr2 / 2.0, nX2[0] )
    X22        = np.linspace( xL2[1] + dt2 / 2.0, xU2[1] - dt2 / 2.0, nX2[1] )
    theta1, r1 = np.meshgrid( X21, X11 )
    theta2, r2 = np.meshgrid( 2.0 * np.pi - X22, X12 )

    nSS1 = FileArray1.shape[0]
    nSS2 = FileArray2.shape[0]

    uK1 = np.empty( (nSS1,nX1[0]), np.float64 )
    uK2 = np.empty( (nSS2,nX2[0]), np.float64 )

    for iSS in range( nSS1 ):

        uK1[iSS] = ComputeAngleAverage( nX1, X22, Data1[iSS], dt1 )
        uK2[iSS] = ComputeAngleAverage( nX2, X22, Data2[iSS], dt2 )

#        den1 = np.copy( uK1[iSS] )
#        den2 = np.copy( uK2[iSS] )
#
#        ind = np.where( np.abs( uK1[iSS] ) < 1.0e-16 )[0]
#        for i in range( ind.shape[0] ):
#            den1[ind[i]] = 1.0e-17
##            uK1[iSS,ind[i]] = 1.0e-17
#
#        ind = np.where( np.abs( uK2[iSS] ) < 1.0e-16 )[0]
#        for i in range( ind.shape[0] ):
#            den2[ind[i]] = 1.0e-17
##            uK2[iSS,ind[i]] = 1.0e-17

        for iX2 in range( nX1[1] ):
            Data1[iSS,:,iX2] = ( Data1[iSS,:,iX2]  )#- uK1[iSS] )# / den1
            Data2[iSS,:,iX2] = ( Data2[iSS,:,iX2]  )#- uK2[iSS] )# / den2

    if( UseCustomTicks1 and UseLogScale1 ):

        ind = np.where( np.abs( Data1 ) < 1.0e-16 )
        nZeros = ind[0].shape[0]
        for i in range( nZeros ):
            Data1[ind[0][i],ind[1][i],ind[2][i]] = 1.0e-17

    if( UseCustomTicks2 and UseLogScale2 ):

        ind = np.where( np.abs( Data2 ) < 1.0e-16 )
        nZeros = ind[0].shape[0]
        for i in range( nZeros ):
            Data2[ind[0][i],ind[1][i],ind[2][i]] = 1.0e-17

    vmin1 = min( +np.inf, Data1.min() )
    vmax1 = max( -np.inf, Data1.max() )
    vmin2 = min( +np.inf, Data2.min() )
    vmax2 = max( -np.inf, Data2.max() )

    if  ( np.abs( vmax1 ) / np.abs( vmin1 ) > 1.0e10 ):
        vmin1 = vmax1 / 1.0e10
    elif( np.abs( vmin1 ) / np.abs( vmax1 ) > 1.0e10 ):
        vmax1 = vmin1 / 1.0e10

    if  ( np.abs( vmax2 ) / np.abs( vmin2 ) > 1.0e10 ):
        vmin2 = vmax2 / 1.0e10
    elif( np.abs( vmin2 ) / np.abs( vmax2 ) > 1.0e10 ):
        vmax2 = vmin2 / 1.0e10

    if( Field1 == Field2 ):

        vmin1 = min( vmin1, vmin2 )
        vmax1 = max( vmax1, vmax2 )
        vmin2 = vmin1
        vmax2 = vmax1

    vmin1 = 1.0e+14
    vmax1 = 1.0e+16
    vmin2 = -1.5e3
    vmax2 = +1.5e3

    if vmin1 < 0.0 and UseLogScale1: UseCustomTicks1 = False
    if vmin2 < 0.0 and UseLogScale2: UseCustomTicks2 = False

    if( UseCustomTicks1 ):

        nTicks1 = 5

        if( UseLogScale1 ):

            ticks1 \
              = np.logspace( np.log10( vmin1 ), np.log10( vmax1 ), nTicks1 )

        else:

            ticks1 = np.linspace( vmin1, vmax1, nTicks1 )

        ticklabels1 = []
        for tick in ticks1:
            ticklabels1.append( '{:.3e}'.format( tick ) )

    if( UseCustomTicks2 ):

        nTicks2 = 5

        if( UseLogScale2 ):

            ticks2 \
              = np.logspace( np.log10( vmin2 ), np.log10( vmax2 ), nTicks2 )

        else:

            ticks2 = np.linspace( vmin2, vmax2, nTicks2 )

        ticklabels2 = []
        for tick in ticks2:
            ticklabels2.append( '{:.3e}'.format( tick ) )

    Norm1 = GetNorm( UseLogScale1, Data1, vmin = vmin1, vmax = vmax1 )
    Norm2 = GetNorm( UseLogScale2, Data2, vmin = vmin2, vmax = vmax2 )

    def f1(t):
        return Data1[t]

    def f2(t):
        return Data2[t]

    # Taken from:
    # https://brushingupscience.com/2016/06/21/
    # matplotlib-animations-the-easy-way/
    im1 = ax.pcolormesh( theta1, r1, f1(0)[:,:], \
                         cmap = cmap1, \
                         norm = Norm1, \
                         shading = 'nearest' )
    im2 = ax.pcolormesh( theta2, r2, f2(0)[:,:], \
                         cmap = cmap2, \
                         norm = Norm2, \
                         shading = 'nearest' )

    # Limits on coordinate axes

    rmax = 360.0
    ax.set_thetamin( 0.0 )
    ax.set_thetamax( 360.0 )
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

    cb1axes = fig.add_axes( [ 0.81,  0.1, 0.03, 0.8 ] )
    cb2axes = fig.add_axes( [ 0.185, 0.1, 0.03, 0.8 ] )
    cbar1 = fig.colorbar( im1, cax = cb1axes )
    cbar2 = fig.colorbar( im2, cax = cb2axes )

    if( UseCustomTicks1 ):

        cbar1.ax.minorticks_off()
        cbar1.set_ticks( ticks1 )
        cbar1.ax.set_yticklabels( ticklabels1 )

    if( UseCustomTicks2 ):

        cbar2.ax.minorticks_off()
        cbar2.set_ticks( ticks2 )
        cbar2.ax.set_yticklabels( ticklabels2 )

    cbar1.set_label( 'Entropy [{:}]'.format( Field1Unit ) )
    cbar2.set_label( 'V2 [{:}]'.format( Field2Unit ) )

    cb2axes.yaxis.set_ticks_position( 'left' )
    cb2axes.yaxis.set_label_position( 'left' )

    def UpdateFrame(t):
        im1.set_array( f1(t)[:,:].flatten() )
        im2.set_array( f2(t)[:,:].flatten() )
        time_text.set_text( 'time = {:d} ms'.format( np.int( Time1[t] ) ) )
        return im1, im2, time_text

    # Call the animator

    print( 'Making movie...' )

    nFrames = min( FileArray1.shape[0], FileArray2.shape[0] )
    fps = max( np.int64( nFrames / MovieRunTime ), 1 )

    anim \
      = animation.FuncAnimation \
          ( fig, UpdateFrame, frames = nFrames, blit = True )

    anim.save( MovieName, fps = fps )
    plt.close()

    os.system( 'rm -rf __pycache__ ' )

Field1 = 'Entropy'
ID1 = 'NR2D_M1.4_Mdot0.3_Rs180_PA1.00e-04_nX640x064'

Field2 = 'PF_V2'
ID2 = 'NR2D_M1.4_Mdot0.3_Rs180_PA1.00e-04_nX640x064'

MovieName \
  = 'mov.{:}_{:}_{:}_{:}.mp4'.format \
      ( ID1, ID2, Field1, Field2 )
MakeMovie2D( ID1, ID2, Field1, Field2, MovieName )
