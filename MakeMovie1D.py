#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from sys import exit
import os

from MakeDataFile import MakeDataFile

def MakeMovie1D( SSi = -1, SSf = -1 ):

    # === User input ===

    nLines = 1

    Field = [ 'AF_P' ]

    global ID
    ID = [ 'NR1D_M0.14_Mdot0.03_Rs180_PA0.00e-00_nX640' ]

#    DataDirectory \
#      = '/home/dunhamsj/AccretionShockData/'
    DataDirectory \
      = '/lump/data/AccretionShockStudy/'

    labels = ID

    yLabel = Field[0]

    UseLogScale = True

    MovieRunTime = 10.0 # [s]

    SaveFileAs = 'mov.{:}_{:}.mp4'.format( ID[0], Field[0] )

    # === End of user input ===

    global Data
    Data = {}

    global Time
    Time = {}

    fig, ax = plt.subplots()

    global ylim
    ylim = [ +np.inf, -np.inf ]

    if( UseLogScale ): ax.set_yscale( 'log' )

    def GetData( i, SSi, SSf ):

        global ID

        iID    = ID[i]
        iField = Field[i]

        PlotFileBaseName = iID + '.plt'

        iDataDirectory = DataDirectory + iID

        DataFileName = '.' + iID + '_' + iField + '.dat'

        ID[i] += '_' + iField

        xL, xU, nX, FileArray \
          = MakeDataFile( iField, iDataDirectory, \
                          DataFileName, PlotFileBaseName, \
                          SSi = SSi, SSf = SSf )

        if( nX[1] > 1 or nX[2] > 1 ):

            print( 'MakeMovie1D.py incompatible with multi-dimensional data.' )
            exit( 'Exiting...' )

        f = open( DataFileName )
        header = f.readline()[16:-2]
        FieldUnit = f.readline()[9:-1]
        DataShape = tuple( [ np.int64( dim ) for dim in header.split( ',' ) ] )

        iTime = list( [ np.float64( t ) \
                        for t in f.readline()[12:-2].split(' ') ] )
        iTime = np.array( iTime )

        iData = np.loadtxt( DataFileName, dtype = np.float64 ).reshape \
                 ( DataShape )

        global ylim
        ylim = [ min( ylim[0], iData.min() ), max( ylim[1], iData.max() ) ]
        #ylim = [ 1.0e14, 1.0e16 ]

        global Data
        Data[ID[i]] = iData

        global Time
        Time[ID[i]] = iTime

        return xL, xU, nX, FileArray

    if SSi < 0: SSi = 0
    if SSf < 0: SSf = 100

    for i in range( nLines ):

        xL, xU, nX, FileArray = GetData( i, SSi, SSf )

        if( i == 0 ):

            xlim  = [ xL[0], xU[0] ]
            Width = xlim[1] - xlim[0]
            dr    = Width / np.float64( nX[0] )
            r     = np.linspace( xlim[0] + dr / 2.0, xlim[1] - dr / 2.0, nX[0] )

            nFiles = FileArray.shape[0]

    ax.set_xlim( xlim )
    ax.set_ylim( ylim )

    ax.set_xlabel( 'Radial Coordinate [km]' )
    ax.set_ylabel( yLabel )

    ax.axvline( 180.0 )

    Height    = ylim[1] - ylim[0]
    time_text = plt.text( xlim[0] + 0.5 * Width, \
                          ylim[0] + 0.7 * Height, '' )

    lines = np.empty( (nLines), object )

    for i in range( nLines ):
        lines[i], = ax.plot( [], [], '.', label = labels[i], \
                             markersize = 1.0, linewidth = 1 )

    def InitializeFrame():
        ret = []
        for i in range( nLines ):
            lines[i].set_data([],[])
            ret.append( lines[i] )
        time_text.set_text('')
        ret.append( time_text )
        ret = tuple( ret )
        return ret

    def UpdateFrame(t):
        global Data
        global Time
        ret = []
        for i in range( nLines ):
            iData = Data[ID[i]]
            y = iData[t]
            lines[i].set_data( r, y )
            ret.append( lines[i] )
            del y
        iTime = Time[ID[i]]
        time_text.set_text( 'time = {:.3e} ms'.format( iTime[t] ) )
        ret.append( time_text )
        ret = ( ret )
        return ret

#    ax.legend( loc = 2 )

    fps = np.float64( nFiles ) / MovieRunTime

    anim = animation.FuncAnimation( fig, UpdateFrame, \
                                    init_func = InitializeFrame, \
                                    frames = nFiles, \
                                    blit = True )

    anim.save( SaveFileAs, fps = fps, dpi = 300 )

    os.system( 'rm -rf __pycache__ ' )

    return

MakeMovie1D( SSi = 0, SSf = 1999 )
