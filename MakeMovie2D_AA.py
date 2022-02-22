#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.integrate import simps
from sys import exit
import os

from MakeDataFile import MakeDataFile

class MakeMovie2D_AA:

    def __init__( self, DataDirectory, ID, Field, UseLogScale = False, \
                  SSi = -1, SSf = -1, MovieRunTime = 30.0, \
                  SaveMovieAs = 'mov.mp4', suffix = '' ):

        self.DataDirectory = DataDirectory
        self.ID = ID
        self.Field = Field
        self.UseLogScale = UseLogScale
        self.SSi = SSi
        self.SSf = SSf
        self.MovieRunTime = MovieRunTime
        self.SaveMovieAs = SaveMovieAs
        self.suffix = suffix

        return

    def ComputeAngleAverage( self, DataIn ):
        return 1.0 / 2.0 * simps( DataIn * np.sin( self.X2 ), x = self.X2 )

    def GetData( self ):

        self.DataDirectory += self.ID + self.suffix

        PlotFileBaseName = self.ID + '.plt'

        DataFileName = '.' + self.ID + '_' + self.Field + '.dat'

        self.ID += '_' + Field

        self.xL, self.xU, self.nX, self.FileArray \
          = MakeDataFile( Field, self.DataDirectory, \

                          DataFileName, PlotFileBaseName, \
                          SSi = self.SSi, SSf = self.SSf )

        f = open( DataFileName )
        header = f.readline()[16:-2]
        FieldUnit = f.readline()[9:-1]
        DataShape = tuple( [ np.int64( dim ) for dim in header.split( ',' ) ] )

        Time = list( [ np.float64( t ) \
                       for t in f.readline()[12:-2].split(' ') ] )
        self.Time = np.array( Time )

        Data = np.loadtxt( DataFileName, dtype = np.float64 ).reshape \
                 ( DataShape )

        self.dX = ( self.xU - self.xL ) / np.float64( self.nX )

        self.X1 = np.linspace( self.xL[0] + 0.5 * self.dX[0], \
                               self.xU[0] - 0.5 * self.dX[0], self.nX[0] )
        self.X2 = np.linspace( self.xL[1] + 0.5 * self.dX[1], \
                               self.xU[1] - 0.5 * self.dX[1], self.nX[1] )

        self.DataAA = np.empty( (DataShape[0],DataShape[1]), np.float64 )
        for i in range( self.DataAA.shape[0] ):
            for j in range( self.DataAA.shape[1] ):
                self.DataAA[i,j] = self.ComputeAngleAverage( Data[i,j] )

        DataFileName = '.GR1D_M0.14_Mdot0.03_Rs180_PA0.00e-00_nX640_{:}.dat'.format( self.Field )

        f = open( DataFileName )
        header = f.readline()[16:-2]
        FieldUnit = f.readline()[9:-1]
        DataShape = tuple( [ np.int64( dim ) for dim in header.split( ',' ) ] )

        self.Data1D = np.loadtxt( DataFileName, dtype = np.float64 ).reshape( DataShape )

        self.xlim = np.array( [ self.xL[0], self.xU[0] ], np.float64 )
        self.ylim = np.array( [ self.DataAA.min(), self.DataAA.max() ], np.float64 )

        return

    def DefineLines( self, ax ):

        self.time_text \
          = ax.text( 0.5, 0.7, '', transform = ax.transAxes )

        self.lines1, \
          = ax.plot( [], [], 'r.', label = '2D (AA)', markersize = 1.0 )
        self.lines2, \
          = ax.plot( [], [], 'b.', label = '1D', markersize = 1.0 )

        ax.set_xlim( self.xlim )
        ax.set_ylim( self.ylim )
        ax.set_ylabel( self.Field )

        ax.legend()

        if self.UseLogScale: ax.set_yscale( 'log' )

        return

    def InitializeFrame( self ):
        self.lines1.set_data([],[])
        self.lines2.set_data([],[])
        self.time_text.set_text('')
        return ( self.lines1, self.lines2, self.time_text )

    def UpdateFrame( self, t ):
        self.lines1.set_data( self.X1, self.DataAA[t] )
        self.lines2.set_data( self.X1, self.Data1D[t] )
        self.time_text.set_text( 'time = {:.3e} ms'.format( self.Time[t] ) )
        return ( self.lines1, self.lines2, self.time_text )

    def MakeMovie( self, fig ):

        nFiles = self.FileArray.shape[0]

        fps = np.float64( nFiles ) / self.MovieRunTime

        anim = animation.FuncAnimation( fig, self.UpdateFrame, \
                                        init_func = self.InitializeFrame, \
                                        frames = nFiles, \
                                        blit = True )

        anim.save( self.SaveMovieAs, fps = fps, dpi = 300 )

        os.system( 'rm -rf __pycache__ ' )

        return

if __name__ == '__main__':

    #DataDirectory = '/home/dunhamsj/AccretionShockData/'
    DataDirectory = '/lump/data/AccretionShockStudy/'
    ID = 'GR2D_M0.14_Mdot0.03_Rs180_PA0.00e-00_nX640x064'
    Field = 'PF_V1'
    UseLogScale = False
    SSi = 0
    SSf = 100
    MovieRunTime = 10.0 # seconds
    SaveMovieAs = 'mov.{:}_{:}.mp4'.format( ID, Field )
    suffix = '_new'

    Movie \
      = MakeMovie2D_AA \
          ( DataDirectory, ID, Field, UseLogScale, SSi, SSf, \
            MovieRunTime, SaveMovieAs, suffix )

    Movie.GetData()

    fig, ax = plt.subplots()

    Movie.DefineLines( ax )

    Movie.MakeMovie( fig )
