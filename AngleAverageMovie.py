#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import os

from MakeDataFile import MakeDataFile
from UtilitiesModule import GetNorm

def ComputeAngleAverage( nX, theta, Data, dtheta ):

    uK = np.empty( (nX[0]), np.float64 )

    for iX1 in range( uK.shape[0] ):

        uK[iX1] = 0.5 * np.sum( Data[iX1] * np.sin( theta ) ) * dtheta

    return uK

def MakeAngleAverageMovie( ID, Field ):

    ############# User input #############

    Root = '/home/dunhamsj/AccretionShockData'
    DataDirectory = Root + '/{:}/'.format( ID )

    UseLogScale = True # Use log scale for field?

    SaveFileAs = 'mov.{:}_{:}_LineOut.mp4'.format( ID, Field )

    MovieRunTime = 10.0 # seconds

    ############# End of user input #############

    PlotFileBaseName = ID + '.plt'

    ID = '{:}_{:}'.format( ID, Field )

    DataFileName = '.{:}.dat'.format( ID )

    xL, xU, nX, FileArray \
      = MakeDataFile( Field, DataDirectory, DataFileName, \
                      PlotFileBaseName, SSi = 0, SSf = 300 )

    dX1 = ( xU[0] - xL[0] ) / np.float64( nX[0] )
    dX2 = ( xU[1] - xL[1] ) / np.float64( nX[1] )

    X1  = np.linspace( xL[0] + dX1 / 2.0, xU[0] - dX1 / 2.0, nX[0] )
    X2  = np.linspace( xL[1] + dX2 / 2.0, xU[1] - dX2 / 2.0, nX[1] )

    f = open( DataFileName )
    header = f.readline()[16:-2]
    FieldUnit = f.readline()[9:-1]
    DataShape = tuple( [ np.int64( dim ) for dim in header.split( ',' ) ] )

    Time = list( [ np.float64( t ) for t in f.readline()[12:-2].split(' ') ] )
    Time = np.array( Time )

    Data = np.loadtxt( DataFileName, dtype = np.float64 ).reshape( DataShape )

    nSS = FileArray.shape[0]

    uK = np.empty( (nSS,nX[0]), np.float64 )

    for iSS in range( nSS ):

        uK[iSS] = ComputeAngleAverage( nX, X2, Data[iSS], dX2 )

    # Plotting

    def f(t):
        return uK[t]

    fig = plt.figure( figsize = (16,9) )
    ax  = fig.add_subplot( 111 )
    ax.grid()

    #if UseLogScale:

    #  if uK.min() < 0.0: ax.set_yscale( 'symlog' )
    #  else:              ax.set_yscale( 'log' )

    xlim = [ xL[0], xU[0] ]
    ylim = [ uK.min(), uK.max() ]

    time_text = plt.text( 0.5, 0.7, '', transform = ax.transAxes )

    line, = ax.plot( [], [], 'k-' )

    def UpdateFrame(t):
        line.set_data( X1, f(t) )
        time_text.set_text( 'time = {:.3e} ms'.format( Time[t] ) )
        return line, time_text

    ax.set_xlim( xlim[0], xlim[1] )
    ax.set_xlabel( 'Radial Coordinate [km]' )

    ax.set_ylim( ylim[0], ylim[1] )
    ax.set_ylabel( '<{:}> [{:}]'.format( Field, FieldUnit ) )

    fps = max( np.int64( nSS / MovieRunTime ), 1 )

    anim = animation.FuncAnimation( fig, UpdateFrame, \
                                    frames = nSS, \
                                    blit = True )

    anim.save( SaveFileAs, fps = fps )
    plt.close()

    os.system( 'rm -rf __pycache__ ' )

ID    = 'NR2D_M1.4_Mdot0.3_Rs180_PA1.00e-04_nX640x064'
Field = 'PF_V2'

MakeAngleAverageMovie( ID, Field )
