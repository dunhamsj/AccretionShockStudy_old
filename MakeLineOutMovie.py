#!/usr/bin/env python3

import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import os

from UtilitiesModule import GetFileArray, ChoosePlotFile

yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

class MakeLineOutMovie:

    def __init__( self, SSi = -1, SSf = -1, nSS = -1 ):

        self.SSi = SSi
        self.SSf = SSf
        self.nSS = nSS

        return

    def GetData( self, X1, PlotFileDirectory, ID, nX ):

        PlotFileBaseName = ID + '.plt'

        PlotFileArray = GetFileArray( PlotFileDirectory, PlotFileBaseName )

        if self.SSi < 0: self.SSi = 0
        if self.SSf < 0: self.SSf = PlotFileArray.shape[0] - 1
        if self.nSS < 0: self.nSS = PlotFileArray.shape[0]

        Data = np.empty( (self.nSS,X1.shape[0],nX[1]), np.float64 )
        Time = np.empty( (self.nSS), np.float64 )

        for i in range( self.nSS ):

            iSS = self.SSi \
                    + np.int64( ( self.SSf - self.SSi + 1 ) / self.nSS ) * i

            PlotFile \
              = ChoosePlotFile( PlotFileDirectory, PlotFileBaseName, \
                                argv = [ 'a', PlotFileArray[iSS] ] )

            ds = yt.load( '{:}'.format( PlotFileDirectory + PlotFile ) )

            Time[i] = ds.current_time.to_ndarray()

            for iX1 in range( X1.shape[0] ):

                oray = ds.ortho_ray( axis = 1, coords = (X1[iX1],0) )

                Data[i,iX1] = oray['AF_P']

        return Data, Time

if __name__ == '__main__':

    SaveFileAs = 'mov.LineOut.mp4'
    MovieRunTime = 10.0 # s

    UseLogScale = False
    ID = 'GR2D_M2.0_Mdot0.3_Rs150'

    PlotFileDirectory = '/lump/data/AccretionShockStudy/' + ID + '/'

    X1 = np.linspace( 41, 149, 10 )
    X1 = np.array( [ 41.0, 100.0, 149.0 ], np.float64 )

    nX = np.array( [ 640, 64, 1 ], np.int64 )

    xL = np.array( [ 40.0, 0.0, 0.0 ], np.float64 )
    xU = np.array( [ 360.0, np.pi, 2.0 * np.pi ], np.float64 )

    MM = MakeLineOutMovie( nSS = -1 )

    Data, Time = MM.GetData( X1, PlotFileDirectory, ID, nX )

    Data = ( Data - Data[0] ) / Data[0]

    nSS = MM.nSS

    dX2 = ( xU[1] - xL[0] ) / np.float64( nX[1] )

    X2 = np.linspace( xL[1] + 0.5 * dX2, xU[1] - 0.5 * dX2, nX[1] )

    fig, ax = plt.subplots( figsize = (10,6) )

    ax.set_xlim( xL[1], xU[1] )
    ax.set_ylim( Data.min(), Data.max() )
    ax.set_xlabel( r'$\theta$' )
    ax.set_ylabel( r'$\left(P\left(t\right)-P\left(0\right)\right)/P\left(0\right)$' )

    if UseLogScale: ax.set_yscale( 'log' )

    time_text = plt.text( 0.5, 0.9, '', transform = ax.transAxes )

    lines = np.empty( X1.shape[0], object )

    for i in range( lines.shape[0] ):

        lines[i], = ax.plot( [], [], label = 'r = {:.1f} km'.format( X1[i] ) )

    def InitializeFrame():

        ret = []

        for i in range( lines.shape[0] ):
            lines[i].set_data([],[])
            ret.append( lines[i] )

        time_text.set_text('')
        ret.append( time_text )
        ret = tuple( ret )

        return ret

    def UpdateFrame(t):

        ret = []

        for i in range( lines.shape[0] ):
            y = Data[t,i]
            lines[i].set_data( X2, y )
            ret.append( lines[i] )

        time_text.set_text( 'time = {:d} ms'.format( np.int64( Time[t] ) ) )
        ret.append( time_text )

        ret = tuple( ret )
        return ret

    ax.legend( loc = 3 )

    fps = max( 1, np.int64( np.float64( nSS ) / MovieRunTime ) )

    anim = animation.FuncAnimation( fig, UpdateFrame, \
                                    init_func = InitializeFrame, \
                                    frames = nSS, \
                                    blit = True )

    anim.save( SaveFileAs, fps = fps, dpi = 300 )

    os.system( 'rm -rf __pycache__ ' )

