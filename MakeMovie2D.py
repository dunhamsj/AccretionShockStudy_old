#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import subprocess
import os

from MakeDataFile import MakeDataFile
from UtilitiesModule import GetNorm

def MakeMovie2D( ID, Field ):

    """
    Generate a movie from data files created in MakeDataFile.py.
    """

    # --- Get user's HOME directory ---
    HOME = subprocess.check_output( ["echo $HOME"], shell = True)
    HOME = HOME[:-1].decode( "utf-8" ) + '/'

    ############# User input #############

    DataDirectory \
      = HOME + 'Research/Data/AccretionShockParameterStudy/{:}/'.format( ID )

    cmap = 'RdBu' # Color scheme for movie

    UseLogScale = False # Use log scale for field?

    MovieRunTime = 10.0 # seconds

    UseCustomTicks = False # Define limits for colorbar

    zAxisVertical = False # Orient z-axis

    ############# End of user input #############

    PlotFileBaseName = ID + '.plt'

    ID = '{:}_{:}'.format( ID, Field )

    DataFileName = '{:}.dat'.format( ID )
    MovieName    = 'mov.{:}.mp4'.format( ID )

    xL, xU, nX, FileArray \
      = MakeDataFile( Field, DataDirectory, DataFileName, PlotFileBaseName )

    if( not DataFileName[0] == '.' ): DataFileName = '.' + DataFileName

    dr = ( xU[0] - xL[0] ) / np.float64( nX[0] )
    dt = ( xU[1] - xL[1] ) / np.float64( nX[1] )

    f = open( DataFileName  )
    header = f.readline()[16:-2]
    FieldUnit = f.readline()[9:-1]
    DataShape = tuple( [ np.int64( dim ) for dim in header.split( ',' ) ] )

    Time = list( [ np.float64( t ) for t in f.readline()[12:-2].split(' ') ] )
    Time = np.array( Time )

    Data = np.loadtxt( DataFileName, dtype = np.float64 ).reshape( DataShape )

    fig = plt.figure( figsize = (16,9) )
    ax  = fig.add_subplot( 111, polar = True )
    X1  = np.linspace( xL[0] + dr / 2.0, xU[0] - dr / 2.0, nX[0] )
    X2  = np.linspace( xL[1] + dt / 2.0, xU[1] - dt / 2.0, nX[1] )
    theta, r = np.meshgrid( X2, X1 )

    if( UseCustomTicks and UseLogScale ):

        ind = np.where( np.abs( Data ) < 1.0e-16 )
        nZeros = ind[0].shape[0]
        for i in range( nZeros ):
            Data[ind[0][i],ind[1][i],ind[2][i]] = 1.0e-17

    vmax = max( -np.inf, Data.max() )
    vmin = min( +np.inf, Data.min() )
    vmin = -1.0e-6
    vmax = +1.0e-6

    if  ( np.abs( vmax ) / np.abs( vmin ) > 1.0e10 ):
        vmin = vmax / 1.0e10
    elif( np.abs( vmin ) / np.abs( vmax ) > 1.0e10 ):
        vmax = vmin / 1.0e10

    if( UseCustomTicks ):

        nTicks = 10

#        vmin = 0.01

        if( UseLogScale ):

            if  ( vmax < 0.0 ):
                ticks \
                  = np.logspace( -np.log10(-vmin), -np.log10(-vmax), nTicks )
            elif( vmin < 0.0 ):
                ticks \
                  = np.logspace( -np.log10(-vmin), +np.log10(+vmax), nTicks )
            else:
                ticks \
                  = np.logspace( +np.log10(+vmin), +np.log10(+vmax), nTicks )

        else:

            ticks = np.linspace( vmin, vmax, nTicks )

        ticklabels = []
        for tick in ticks:
            ticklabels.append( '{:.3e}'.format( tick ) )

    Norm = GetNorm( UseLogScale, Data, vmin = vmin, vmax = vmax )

    def f(t):
        return Data[t]

    # Taken from:
    # https://brushingupscience.com/2016/06/21/
    # matplotlib-animations-the-easy-way/
    im = ax.pcolormesh( theta, r, f(0)[:,:], \
                        cmap = cmap, \
                        norm = Norm, \
                        shading = 'nearest' )

    # Limits on coordinate axes

    ax.set_thetamin( 180.0/np.pi * xL[1]  )
    ax.set_thetamax( 180.0/np.pi * xU[1] )
#    ax.set_rlim( 0.0, 200.0 )
    ax.set_theta_direction( -1 )

    if( zAxisVertical ):
        ax.set_theta_zero_location( 'N' ) # z-axis vertical
        time_text = plt.text( 0.5 * np.pi / 2, xU[0] * ( 1.0 + 0.3 ), '' )
    else:
        ax.set_theta_zero_location( 'W' ) # z-axis horizontal
        time_text = ax.text( 0.475, 0.7, '', transform = ax.transAxes )

    ax.set_position( [0.1,-0.45,0.7,2] )

    cax = fig.add_axes( [0.85,0.1,0.03,0.8] )
    cbar = fig.colorbar( im, cax = cax )

    if( UseCustomTicks ):

        cbar.ax.minorticks_off()
        cbar.set_ticks( ticks )
        cbar.ax.set_yticklabels( ticklabels )

    cbar.set_label( Field + ' ' + FieldUnit )

    TimeUnit = 'ms'

    def UpdateFrame(t):

        im.set_array( f(t)[:,:].flatten() )
        time_text.set_text( 'Time = {:.6e} {:}'.format \
                              ( Time[t], TimeUnit ) )

        ret = ( im, time_text )

        return ret

    # Call the animator

    print( 'Making movie...' )

    nFrames = FileArray.shape[0]
    fps = nFrames / MovieRunTime

    anim \
      = animation.FuncAnimation \
          ( fig, UpdateFrame, \
            frames = nFrames, blit = True )

    anim.save( MovieName, fps = fps )
    plt.close()

    os.system( 'rm -rf __pycache__ ' )

Field = 'DivV2'
ID = 'GR2D_M1.4_Mdot0.3_Rs180_PA0.000_nX640x064'
MakeMovie2D( ID, Field )
