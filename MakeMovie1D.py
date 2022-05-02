#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from sys import exit
import os

from MakeDataFile_new import MakeDataFile_new

class MakeMovie1D:

    def __init__( self, SSi, SSf, nSS = -1 ):

        self.SSi = SSi
        self.SSf = SSf
        if nSS < 0:
            self.nSS = SSf - SSi + 1
        else:
            self.nSS = nSS

        return


    def GetData( self, DataDirectory, ID, PlotFileBaseName, Field ):

        DataFileDirectory = '.{:}_{:}_MovieData1D'.format( ID, Field )

        xL, xU, nX, FileArray \
          = MakeDataFile_new( Field, DataDirectory + ID + '/', \
                              DataFileDirectory, PlotFileBaseName, \
                              CoordinateSystem = 'cartesian', \
                              SSi = self.SSi, SSf = self.SSf, nSS = self.nSS )
        exit()

#        f = open( DataFileName )
#        header = f.readline()[16:-2]
#        FieldUnit = f.readline()[9:-1]
#        DataShape = tuple( [ np.int64( dim ) for dim in header.split( ',' ) ] )
#
#        Time = list( [ np.float64( t ) \
#                        for t in f.readline()[12:-2].split(' ') ] )
#        Time = np.array( Time )
#
#        Data = np.loadtxt( DataFileName, dtype = np.float64 ).reshape \
#                 ( DataShape )
#
#        xlim  = [ xL[0], xU[0] ]
#        dr    = ( xU[0] - xL[0] ) / np.float64( nX[0] )
#        r     = np.linspace( xlim[0] + 0.5 * dr, xlim[1] - 0.5 * dr, nX[0] )
#
#        self.xlim = np.array( xlim )
#        self.r    = r
#        self.Time = Time
#
#        return Data

if __name__ == '__main__':

    MovieRunTime = 30.0 # [s]
    SaveFileAs   = 'mov.1D.mp4'
    UseLogScale  = True
    SSi          = 0
    SSf          = 1999
    nSS          = 500

    MM1D = MakeMovie1D( SSi = SSi, SSf = SSf, nSS = nSS )

    DataDirectory    = '/lump/data/AccretionShockStudy/'
    ID               = 'GR1D_M2.0_Mdot0.3_Rs150_entropyPert'
    PlotFileBaseName = ID + '.plt'
    MM1D.GetData( DataDirectory, ID, PlotFileBaseName, 'PF_D' )
    exit()
    #PF_D = MM1D.GetData( DataDirectory, ID, PlotFileBaseName, 'PF_D' )

#    PF_D = np.copy( ( PF_D - PF_D[0] ) / PF_D[0] )

    # Plotting

    fig, ax = plt.subplots()

    r    = MM1D.r
    xmin = MM1D.xlim[0]
    xmax = 151.0#MM1D.xlim[1]
    ymin = PF_D.min()
    ymax = 0.01#PF_D.max()

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
