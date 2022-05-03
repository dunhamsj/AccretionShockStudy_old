#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import os

from UtilitiesModule import GetFileArray
from MakeDataFile_new import MakeDataFile_new, ReadHeader

class MakeMovie1D:

    def __init__( self, SSi = -1, SSf = -1, nSS = -1 ):

        self.SSi = SSi
        self.SSf = SSf
        self.nSS = nSS

        return


    def GetData( self, PlotFileDirectory, ID, \
                 PlotFileBaseName, Field, nX1 = 640 ):

        DataFileDirectory = '.{:}_{:}_MovieData1D'.format( ID, Field )

        MakeDataFile_new( Field, PlotFileDirectory + ID + '/', \
                          DataFileDirectory, PlotFileBaseName, \
                          CoordinateSystem = 'cartesian', \
                          SSi = self.SSi, SSf = self.SSf, nSS = self.nSS, \
                          Verbose = True )

        PlotFileArray = GetFileArray( PlotFileDirectory + ID + '/', \
                                      PlotFileBaseName )

        if self.SSi < 0: self.SSi = 0
        if self.SSf < 0: self.SSf = PlotFileArray.shape[0]
        if self.nSS < 0: self.nSS = self.SSf - self.SSi + 1

        Data = np.empty( (self.nSS,nX1), np.float64 )
        Time = np.empty( (self.nSS)    , np.float64 )

        for i in range( self.nSS ):

            iSS = self.SSi + np.int64( ( self.SSf - self.SSi ) / self.nSS ) * i

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
    SaveFileAs   = 'mov.1D.mp4'
    UseLogScale  = False
    nSS          = 500

    MM1D = MakeMovie1D( nSS = nSS )

    PlotFileDirectory = '/lump/data/AccretionShockStudy/'
    ID                = 'GR1D_M2.0_Mdot0.3_Rs150_entropyPert'
    PlotFileBaseName  = ID + '.plt'
    PF_D \
      = MM1D.GetData( PlotFileDirectory, ID, PlotFileBaseName, 'PF_D', \
                      nX1 = 640 )

    PF_D = np.copy( np.diff( PF_D ) )

    # Plotting

    fig, ax = plt.subplots()

    r    = MM1D.X1_C[:-1]
    xmin = MM1D.xlim[0]
    xmax = 151.0#MM1D.xlim[1]
    ymin = PF_D.min()
    ymax = PF_D.max()

    if UseLogScale: ax.set_yscale( 'log' )

    ax.set_xlabel( 'Radial Coordinate [km]' )
    ax.set_ylabel( 'PF_D' )
    ax.set_xlim( xmin, xmax )
    ax.set_ylim( ymin, ymax )

    line, = ax.plot( [], [], 'r-' )
    time_text = plt.text( 0.5, 0.7, '', transform = ax.transAxes )

    def InitializeFrame():
        line.set_data([],[])
        time_text.set_text('')
        return ( line, time_text )

    def UpdateFrame(t):
        line.set_data( r, PF_D[t] )
        time_text.set_text( 'time = {:.3e} ms'.format( MM1D.Time[t] ) )
        return ( line, time_text )

    fps = max( 1, np.float64( nSS ) / MovieRunTime )

    anim = animation.FuncAnimation( fig, UpdateFrame, \
                                    init_func = InitializeFrame, \
                                    frames = nSS, \
                                    blit = True )

    anim.save( SaveFileAs, fps = fps, dpi = 300 )

    os.system( 'rm -rf __pycache__ ' )
