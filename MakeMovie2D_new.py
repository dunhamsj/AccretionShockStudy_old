#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import os

from MakeDataFile_new import MakeDataFile, ReadHeader

from UtilitiesModule import GetNorm, GetFileArray

class MakeMovie2D:

    def __init__( self, SSi = -1, SSf = -1, nSS = -1 ):

        self.SSi = SSi
        self.SSf = SSf
        self.nSS = nSS

        return

    def ComputeAngleAverage( self, nX, theta, Data, dtheta ):

        uK = np.empty( (nX[0]), np.float64 )

        for iX1 in range( uK.shape[0] ):

            uK[iX1] = 0.5 * np.sum( Data[iX1] * np.sin( theta ) ) * dtheta[0]

        return uK

    def GetData( self, PlotFileDirectory, ID, \
                 PlotFileBaseName, Field, nX1 = 640, nX2 = 64 ):

        DataFileDirectory = '.{:}_{:}_MovieData2D'.format( ID, Field )

        MakeDataFile( Field, PlotFileDirectory + ID + '/', \
                      DataFileDirectory, PlotFileBaseName, \
                      CoordinateSystem = 'spherical', \
                      SSi = self.SSi, SSf = self.SSf, nSS = self.nSS, \
                      Verbose = True )

        PlotFileArray = GetFileArray( PlotFileDirectory + ID + '/', \
                                      PlotFileBaseName )

        if self.SSi < 0: self.SSi = 0
        if self.SSf < 0: self.SSf = PlotFileArray.shape[0] - 1
        if self.nSS < 0: self.nSS = PlotFileArray.shape[0]

        Data = np.empty( (self.nSS,nX1,nX2), np.float64 )
        Data_K = np.empty( (self.nSS,nX1), np.float64 )
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

            Data_K[i] \
              = self.ComputeAngleAverage( DataShape, X2_C, Data[i], dX2 )

        self.xL = np.array( [ X1_C[0 ] - 0.5 * dX1[0 ], \
                              X2_C[0 ] - 0.5 * dX2[0 ] ], np.float64 )
        self.xU = np.array( [ X1_C[-1] + 0.5 * dX1[-1], \
                              X2_C[-1] + 0.5 * dX2[-1] ], np.float64 )

        self.DataUnits = DataUnits
        self.theta, self.r = np.meshgrid( X2_C, X1_C )

        self.rlim     = np.array( [ self.xL[0], self.xU[0] ], np.float64 )
        self.thetalim = np.array( [ self.xL[1], self.xU[1] ], np.float64 )

        self.Time = Time

        return Data, Data_K

if __name__ == '__main__':

    MovieRunTime  = 30.0 # [s]
    SaveFileAs    = 'mov.2D.mp4'
    UseLogScale   = True
    cmap          = 'viridis'
    zAxisVertical = False
    nSS           = -1
    Field         = 'AF_P'

    MM2D = MakeMovie2D( nSS = nSS )

    PlotFileDirectory = '/lump/data/accretionShockStudy/'
    ID                = 'GR2D_M2.0_Mdot0.3_Rs150_ColdStart'
    PlotFileBaseName  = ID + '.plt_'
    Data, Data_K \
      = MM2D.GetData( PlotFileDirectory, ID, PlotFileBaseName, Field, \
                      nX1 = 640, nX2 = 64 )

    if nSS < 0: nSS = MM2D.nSS

    Time      = MM2D.Time
    DataUnits = MM2D.DataUnits
    xL        = MM2D.xL
    xU        = MM2D.xU
    theta     = MM2D.theta
    r         = MM2D.r

#    for t in range( nSS ):
#
#        for iX2 in range( 64 ):
#
#            Data[t,:,iX2] = ( Data[t,:,iX2] - Data_K[t] ) / Data_K[t]

    vmax = Data[0].max()
    vmin = Data[0].min()

    Norm = GetNorm( UseLogScale, Data[0], vmin = vmin, vmax = vmax )#, \
#                    linthresh = 1.0e-8 )

    def f(t):
        return Data[t]

    fig = plt.figure( figsize = (16,9) )
    ax  = fig.add_subplot( 111, polar = True )

    ax.grid( False )

    # Taken from:
    # https://brushingupscience.com/2016/06/21/
    # matplotlib-animations-the-easy-way/
    im = ax.pcolormesh( theta, r, f(0)[:,:], \
                        cmap = cmap, \
                        norm = Norm, \
                        shading = 'nearest' )

    # Limits on coordinate axes

    ax.set_thetamin( 180.0/np.pi * xL[1] )
    ax.set_thetamax( 180.0/np.pi * xU[1] )
   # ax.set_rlim( 80.0, 100.0 )
   # ax.set_ylim( [80.0,100.0] )
    ax.set_rorigin(0.0)
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

    cbar.set_label( Field + ' ' + DataUnits )

    TimeUnit = 'ms'

    def UpdateFrame(t):

        im.set_array( f(t)[:,:].flatten() )
        time_text.set_text( 'Time = {:.6e} {:}'.format \
                              ( Time[t], TimeUnit ) )

        return ( im, time_text )

    # Call the animator

    print( 'Making movie...' )

    fps = max( 1, np.int64( nSS / MovieRunTime ) )

    anim \
      = animation.FuncAnimation( fig, UpdateFrame, frames = nSS, blit = True )

    anim.save( SaveFileAs, fps = fps, dpi = 300 )

    plt.close()

    os.system( 'rm -rf __pycache__ ' )
