#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from sys import exit
import os

from MakeDataFile import MakeDataFile

def MakeMovie1D():

    # === User input ===

    Field = 'DivV2'

    ID = 'GR2D_M1.4_Mdot0.3_Rs180_PA0.000_nX320x064'

    DataDirectory \
      = '/Users/dunhamsj/Research/Data/AccretionShockParameterStudy/'

    MovieRunTime = 10.0 # s

    SaveFileAs = 'mov.{:}_{:}.mp4'.format( ID, Field )

    RadiiBelowShock = np.array( [ 100.0 ], np.float64 ) # km
    RadiiAboveShock = np.array( [ 300.0 ], np.float64 ) # km

    # === End of user input ===

    fig, axs = plt.subplots( 2, 1 )

    DataFileName = '.' + ID + '_' + Field + '.dat'

    PlotFileBaseName = ID + '.plt'

    DataDirectory += ID + '/'

    xL, xU, nX, FileArray \
      = MakeDataFile( Field, DataDirectory, \
                      DataFileName, PlotFileBaseName )

    f = open( DataFileName )
    header = f.readline()[16:-2]
    FieldUnit = f.readline()[9:-1]
    DataShape = tuple( [ np.int64( dim ) for dim in header.split( ',' ) ] )

    Time = list( [ np.float64( t ) \
                   for t in f.readline()[12:-2].split(' ') ] )
    Time = np.array( Time )

    Data = np.loadtxt( DataFileName, dtype = np.float64 ).reshape \
            ( DataShape )

    dX = ( xU - xL ) / np.float64( nX )

    r     = np.linspace( xL[0] + dX[0] / 2.0, xU[0] - dX[0] / 2.0, nX[0] )
    theta = np.linspace( xL[1] + dX[1] / 2.0, xU[1] - dX[1] / 2.0, nX[1] )

    xlim  = [ xL[1], xU[1] ]
    Width = xlim[1] - xlim[0]

    indT = np.linspace( 0, 300, 301, dtype = np.int64 )

    ind0 = np.empty( (RadiiBelowShock.shape[0]), np.int64 )
    ind1 = np.empty( (RadiiAboveShock.shape[0]), np.int64 )

    for i in range( ind0.shape[0] ):
        ind0[i] \
          = np.where(  ( r > RadiiBelowShock[i] - 0.75 * dX[0] ) \
                     & ( r < RadiiBelowShock[i] + 0.75 * dX[0] ) )[0][0]

    for i in range( ind1.shape[0] ):
        ind1[i] \
          = np.where(   ( r > RadiiAboveShock[i] - 0.75 * dX[0] ) \
                      & ( r < RadiiAboveShock[i] + 0.75 * dX[0] ) )[0][0]

    ylim0 = [ Data[indT,ind0].min(), Data[indT,ind0].max() ]
    ylim1 = [ Data[indT,ind1].min(), Data[indT,ind1].max() ]
    Height = ylim0[1] - ylim0[0]

    nFiles = indT.shape[0]#FileArray.shape[0]

    axs[0].set_xlim( xlim )
    axs[1].set_xlim( xlim )

    axs[0].set_ylim( ylim0 )
    axs[1].set_ylim( ylim1 )

    axs[0].set_xlabel( r'$\theta$' )
    axs[1].set_xlabel( r'$\theta$' )

    axs[0].set_ylabel( '{:} {:}'.format( Field, FieldUnit ) )
    axs[1].set_ylabel( '{:} {:}'.format( Field, FieldUnit ) )

    time_text = axs[0].text( xlim [0] + 0.6 * Width, \
                             ylim0[0] + 0.9 * Height, '' )

    lines0 = np.empty( ind0.shape[0], object )
    lines1 = np.empty( ind1.shape[0], object )

    for i in range( lines0.shape[0] ):
        lines0[i], = axs[0].plot( [], [], 'k.', markersize = 1.0, \
                          label = 'r = {:.1f} km'.format( r[ind0[i]] ) )

    for i in range( lines1.shape[0] ):
        lines1[i], = axs[1].plot( [], [], 'k.', markersize = 1.0, \
                          label = 'r = {:.1f} km'.format( r[ind1[i]] ) )


    axs[0].legend( loc = 2 )
    axs[1].legend( loc = 2 )

    def InitializeFrame():

        ret = []

        for i in range( lines0.shape[0] ):
            lines0[i].set_data([],[])
            ret.append( lines0[i] )

        for i in range( lines1.shape[0] ):
            lines1[i].set_data([],[])
            ret.append( lines1[i] )

        time_text.set_text('')
        ret.append( time_text )
        ret = tuple( ret )

        return ret

    def UpdateFrame(t):

        ret = []

        for i in range( lines0.shape[0] ):
            y = Data[t,ind0[i]]
            lines0[i].set_data( theta, y )
            ret.append( lines0[i] )

        for i in range( lines1.shape[0] ):
            y = Data[t,ind1[i]]
            lines1[i].set_data( theta, y )
            ret.append( lines1[i] )

        time_text.set_text( 'time = {:d} ms'.format( np.int64( Time[t] ) ) )
        ret.append( time_text )

        ret = tuple( ret )
        return ret

    fps = np.float64( nFiles ) / MovieRunTime

    plt.subplots_adjust( hspace = 0.4 )

    anim = animation.FuncAnimation( fig, UpdateFrame, \
                                    init_func = InitializeFrame, \
                                    frames = nFiles, \
                                    blit = True )

    anim.save( SaveFileAs, fps = fps, dpi = 300 )

    os.system( 'rm -rf __pycache__ ' )

    return

MakeMovie1D()
