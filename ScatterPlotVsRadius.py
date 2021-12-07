#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

class ScatterPlot:

    def __init__( self, Root, ID, nTheta, Field, UseLogScale, \
                  t0, t1, r0, r1, RunTime ):

        print( '\nCreating instance of ScatterPlot class...\n' )

        self.Root        = Root
        self.ID          = ID
        self.nTheta      = nTheta
        self.Field       = Field
        self.UseLogScale = UseLogScale
        self.t0          = t0
        self.t1          = t1
        self.r0          = r0
        self.r1          = r1
        self.RunTime     = RunTime

        self.SaveFileAs \
          = 'mov.{:}_Scatter_{:}.mp4'.format( self.ID, self.Field )

        print( '  Variables:' )
        print( '  ----------' )
        print( '    Root       : {:s}'.format( self.Root ) )
        print( '    ID         : {:s}'.format( self.ID ) )
        print( '    nTheta     : {:d}'.format( self.nTheta ) )
        print( '    Field      : {:s}'.format( self.Field ) )
        print( '    UseLogScale: {:}'.format( self.UseLogScale ) )
        print( '    t0         : {:.3e} ms'.format( self.t0 ) )
        print( '    t1         : {:.3e} ms'.format( self.t1 ) )
        print( '    r0         : {:.3e} km'.format( self.r0 ) )
        print( '    r1         : {:.3e} km'.format( self.r1 ) )
        print( '    RunTime    : {:.3e} s'.format( self.RunTime ) )
        print( '    SaveFileAs : {:s}\n'.format( self.SaveFileAs ) )

        return


    def GetData( self ):

        from MakeDataFile import MakeDataFile

        PlotFileBaseName = self.ID + '.plt'
        DataDirectory    = self.Root + self.ID + '/'
        DataFileName     = '.{:}_{:}.dat'.format( self.ID, self.Field )

        self.xL, self.xU, self.nX, self.FileArray \
          = MakeDataFile( self.Field, DataDirectory, DataFileName, \
                          PlotFileBaseName, UsePhysicalUnits = True, \
                          WriteExtras = False )

        f = open( DataFileName )
        header = f.readline()[16:-2]
        self.FieldUnit = f.readline()[9:-1]
        DataShape = tuple( [ np.int64( dim ) for dim in header.split( ',' ) ] )

        Time = list( [ np.float64( t ) \
                       for t in f.readline()[12:-2].split(' ') ] )
        self.Time = np.array( Time )

        self.AllData = np.loadtxt( DataFileName, dtype = np.float64 ).reshape \
                                   ( DataShape )

        self.dX = ( self.xU - self.xL ) / np.float64( self.nX )

        self.r \
          = np.linspace( self.xL[0] + 0.5 * self.dX[0], \
                         self.xU[0] - 0.5 * self.dX[0], self.nX[0] )

        self.theta \
          = np.linspace( self.xL[1] + 0.5 * self.dX[1], \
                         self.xU[1] - 0.5 * self.dX[1], self.nX[1] )

        return


    def GetPlotIndices( self ):

        self.GetData()

        self.indT  = np.where(   ( self.Time >= self.t0 ) \
                               & ( self.Time <= self.t1 ) )[0]
        self.indX1 = np.where(   ( self.r    >= self.r0 ) \
                               & ( self.r    <= self.r1 ) )[0]

        self.indX2 = np.linspace( 0, self.AllData.shape[2]-1, self.nTheta, \
                                  dtype = np.int64 )

        self.Data = np.empty( (self.indT.shape[0], \
                               self.indX1.shape[0], \
                               self.indX2.shape[0]), np.float64 )

        for i in range( self.indT.shape[0] ):
            for j in range( self.indX1.shape[0] ):
                for k in range( self.indX2.shape[0] ):

                    self.Data[i,j,k] = self.AllData[self.indT[i], \
                                                    self.indX1[j], \
                                                    self.indX2[k]]

        return


    def MakeScatterPlot( self ):

        self.GetPlotIndices()

        xlim  = [ self.r[self.indX1[0]], self.r[self.indX1[-1]] ]
        Width = xlim[1] - xlim[0]

        yMin = +np.inf
        yMax = -np.inf

        for i in range( self.indT.shape[0] ):
            for j in range( self.indX1.shape[0] ):
                yMin = min( yMin, self.Data[i,j].min() )
                yMax = max( yMax, self.Data[i,j].max() )

        if self.UseLogScale and yMin < 0.0:
            yMin, yMax = min( yMin, -yMax ), max( yMax, -yMin )

        yMin *= 1.0
        yMax *= 1.0

        yMin = -1.0e-6#0.005
        yMax = -yMin

        ylim   = [ yMin, yMax ]
        Height = ylim[1] - ylim[0]

        fig, ax = plt.subplots()

        ax.set_xlim( xlim )
        ax.set_ylim( ( ylim[0], ylim[1] ) )

        ax.set_xlabel( 'Radial Coordinate [km]' )
        ax.set_ylabel( self.Field + ' ' + self.FieldUnit )

        lines = np.empty( (self.indX2.shape[0]), object )

        alpha = np.linspace( 0.2, 1.0, self.indX2.shape[0] )

        for i in range( self.indX2.shape[0] ):
            lines[i], \
              = ax.plot( [], [], 'k.', markersize = 2.0, alpha = alpha[i], \
                         label = r'$\theta=$ {:.0f}$\degree$'.format \
                                 ( 180.0 / np.pi * self.theta[self.indX2[i]] ) )

        time_text = ax.text( 0.5, 0.7, '', transform = ax.transAxes )

        def InitializeFrame():

            ret = []

            for i in range( self.indX2.shape[0] ):
                lines[i].set_data([],[])
                ret.append( lines[i] )

            time_text.set_text('')
            ret.append( time_text )

            ret = ( ret )
            return ret

        def UpdateFrame(t):

            ret = []

            for i in range( self.indX2.shape[0] ):
                lines[i].set_data( self.r[self.indX1], self.Data[t,:,i] )

                ret.append( lines[i] )

            time_text.set_text( 'time = {:d} ms'.format \
              ( np.int64( self.Time[self.indT[t]] ) ) )
            ret.append( time_text )

            ret = ( ret )
            return ret

        if self.UseLogScale:
            if yMin < 0.0:
                ax.set_yscale( 'symlog', linthresh = 1.0e-6, base = 10 )
            else:
                ax.set_yscale( 'log' )

        ax.legend(loc=3)

        ax.set_xlabel( 'Radial Coordinate [km]' )
        ax.set_ylabel( '{:} {:}'.format( self.Field, self.FieldUnit ) )

        fps = max( 1, np.int64( np.float64( self.indT.shape[0] ) \
                / self.RunTime ) )

        anim = animation.FuncAnimation( fig, UpdateFrame, \
                                        init_func = InitializeFrame, \
                                        frames = self.indT.shape[0], \
                                        blit = True )

        anim.save( self.SaveFileAs, fps = fps, dpi = 300 )


if __name__ == '__main__':

    Root        = '/Users/dunhamsj/Research/Data/AccretionShockParameterStudy/'
    ID          = 'GR2D_M1.4_Mdot0.3_Rs180_PA0.000_nX640x064'
    nTheta      = 8
    Field       = 'DivV2'
    UseLogScale = False
    t0          = 0.0
    t1          = 600.0
    r0          = 160.0
    r1          = 200.0
    RunTime     = 10.0

    SP = ScatterPlot( Root, ID, nTheta, Field, UseLogScale, \
                      t0, t1, r0, r1, RunTime )

    SP.MakeScatterPlot()

    import os
    os.system( 'rm -rf __pycache__ ' )
