#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import os

from UtilitiesModule import GetFileArray
from MakeDataFile_new import MakeDataFile_new, ReadHeader

class MakeMovie1D:

    def __init__( self, SSi = -1, SSf = -1, nSS = -1, suffix = '' ):

        self.SSi = SSi
        self.SSf = SSf
        self.nSS = nSS
        self.suffix = suffix

        return


    def GetData( self, PlotFileDirectory, ID, \
                 PlotFileBaseName, Field, nX1 = 640 ):

        DataFileDirectory \
          = '.{:}_{:}{:}_MovieData1D'.format( ID, Field, self.suffix )

        PlotFileDirectory += ID + self.suffix + '/'

        MakeDataFile_new( Field, PlotFileDirectory, \
                          DataFileDirectory, PlotFileBaseName, \
                          CoordinateSystem = 'spherical', \
                          SSi = self.SSi, SSf = self.SSf, nSS = self.nSS, \
                          Verbose = True )

        PlotFileArray = GetFileArray( PlotFileDirectory, \
                                      PlotFileBaseName )

        if self.SSi < 0: self.SSi = 0
        if self.SSf < 0: self.SSf = PlotFileArray.shape[0] - 1
        if self.nSS < 0: self.nSS = PlotFileArray.shape[0]

        Data = np.empty( (self.nSS,nX1), np.float64 )
        Time = np.empty( (self.nSS)    , np.float64 )

        for i in range( self.nSS ):

            iSS = self.SSi \
                    + np.int64( ( self.SSf - self.SSi + 1 ) / self.nSS ) * i

            PlotFile = PlotFileArray[iSS]

            DataFile = DataFileDirectory + '/' + PlotFile + '.dat'

            DataShape, DataUnits, t, X1_C, X2_C, X3_C, dX1, dX2, dX3 \
              = ReadHeader( DataFile )

            Time[i] = t
            Data[i] = np.loadtxt( DataFile )

        xL = X1_C[0 ] - 0.5 * dX1[0 ]
        xU = X1_C[-1] + 0.5 * dX1[-1]

        self.X1_C = X1_C

        self.xlim = np.array( [ xL, xU ], np.float64 )

        self.Time = Time

        return Data

if __name__ == '__main__':

    MovieRunTime = 30.0 # [s]
    UseLogScale  = True
    nSS          = -1
    Field        = 'PF_D'
    suffix       = ''

    MM1D = MakeMovie1D( nSS = nSS, suffix = suffix )

    PlotFileDirectory = '/lump/data/AccretionShockStudy/'
    IDs               = np.array( \
      [ 'GR1D_M1.4_Mdot0.3_Rs180' ], str )

    Data = np.empty( IDs.shape[0], object )

    ymin = +np.inf
    ymax = -np.inf

    for i in range( Data.shape[0] ):
        PlotFileBaseName = IDs[i] + '.plt'
        d = MM1D.GetData( PlotFileDirectory, IDs[i], \
                          PlotFileBaseName, Field, nX1 = 640 )
        BG = np.copy( d[-1] )
        Data[i] = np.copy(d)#( d[:-1] - BG ) / BG

        ymin = min( ymin, Data[i].min() )
        ymax = max( ymax, Data[i].max() )

    SaveFileAs   = 'mov.{:}_1D{:}.mp4'.format( IDs[0][5:], suffix )

    nSS = MM1D.nSS# - 1

    # Plotting

    fig, ax = plt.subplots()

    fig.suptitle( IDs[0][5:] )

    r    = MM1D.X1_C
    xmin = MM1D.xlim[0]
    xmax = MM1D.xlim[1]
    ymin = Data.min()
    ymax = Data.max()
    plt.plot( r, Data[0,0] );plt.show();exit()

    if UseLogScale: ax.set_yscale( 'log' )

    ax.set_xlabel( 'Radial Coordinate [km]' )
    ax.set_ylabel( r'$\left(\rho\left(t\right)-\rho\left(0\right)\right)/\rho\left(0\right)$' )
    #ax.set_xlim( xmin, xmax )
    #ax.set_ylim( ymin, ymax )
    ax.grid()

    c = np.array( [ 'r-', 'b-', 'm-' ] )

    lines = np.empty( Data.shape[0], object )
    for i in range( lines.shape[0] ):
        lines[i], = ax.plot( [], [], c[i], label = IDs[i][0:2] )
    time_text = plt.text( 0.2, 0.9, '', transform = ax.transAxes )

    ax.legend( loc = 3 )

    def InitializeFrame():
        ret = []
        for i in range( lines.shape[0] ):
            lines[i].set_data( [], [] )
            ret.append( lines[i] )
        time_text.set_text('')
        ret.append( time_text )
        return tuple( ret )

    def UpdateFrame(t):
        ret = []
        for i in range( lines.shape[0] ):
            lines[i].set_data( r, Data[i][t] )
            ret.append( lines[i] )
        time_text.set_text( 'time = {:.3e} ms'.format( MM1D.Time[t] ) )
        ret.append( time_text )
        return tuple( ret )

    fps = max( 1, np.float64( nSS ) / MovieRunTime )

    anim = animation.FuncAnimation( fig, UpdateFrame, \
                                    init_func = InitializeFrame, \
                                    frames = nSS, \
                                    blit = True )

    anim.save( SaveFileAs, fps = fps, dpi = 300 )

    os.system( 'rm -rf __pycache__ ' )
