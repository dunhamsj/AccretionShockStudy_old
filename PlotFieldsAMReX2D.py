#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
from sys import argv

from UtilitiesModule import GetData as GD
from UtilitiesModule import GetNorm

#Root = '/home/dunhamsj/AccretionShockData/'
Root = '/lump/data/accretionShockStudy/'
#Root = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Euler_Relativistic_IDEAL/'

class PlotFieldsAMReX2D:

    def __init__( self, ID, Field, UseLogScale = True, \
                  suffix = '', cmap = 'viridis' ):

        self.ID = ID
        self.Field = Field

        self.DataDirectory = Root + ID + '{:}/'.format( suffix )
        self.UseLogScale = UseLogScale

        self.cmap = cmap

        return

    def GetData( self ):

        PlotFileBaseName = self.ID + '.plt_'

        self.Data, self.DataUnit, self.r, self.theta, x3, dx1, dx2, dx3, \
          self.xL, self.xU, nX, self.Time \
            = GD( self.DataDirectory, PlotFileBaseName, \
                  self.Field, CoordinateSystem = 'spherical', \
                  UsePhysicalUnits = True, argv = argv, \
                  ReturnMesh = True, ReturnTime = True, Verbose = True )

        self.nX = np.array( [ self.Data.shape[0], self.Data.shape[1], 1 ], \
                            np.int64 )

        self.dX = ( self.xU - self.xL ) / np.float64( self.nX )

        return


    def ComputeAngleAverage( self ):

        self.uK = np.empty( (self.nX[0]), np.float64 )

        theta = self.theta

        for iX1 in range( self.nX[0] ):

            self.uK[iX1] \
              = 0.5 * np.sum( self.Data[iX1] * np.sin( theta ) ) * self.dX[1]

        return

    def PlotWhat( self, ys, yN, iX2, UseLogScale ):

        if UseLogScale:
            plt.semilogy( self.r[:,iX2], \
                          np.abs( self.Data[:,iX2] - ys ) / yN, \
                          label = 'iX2 = {:}'.format( iX2 ) )
        else:
            plt.plot( self.r[:,iX2], \
                          np.abs( self.Data[:,iX2] - ys ) / yN, \
                          label = 'iX2 = {:}'.format( iX2 ) )

        return

    def PlotData( self ):

        UseLogScale = False
#        ys = 0.0
#        yN = 1.0
#        iX2 = np.array( [ 0, 16, 32, 48, 63 ], np.int64 )
#        for i in range( iX2.shape[0] ):
#            self.PlotWhat( ys, yN, iX2[i], UseLogScale )
#        plt.legend()
#        plt.suptitle( self.Field )
#        plt.show()
#        exit()

        fig = plt.figure( figsize = (16,9) )
        ax  = fig.add_subplot( 111, polar = True )

        fig.suptitle( self.ID + '\nTime = ' \
                      + str( np.int64( self.Time ) ) + ' ms' )

        ax.set_thetamin( 0.0  )
        ax.set_thetamax( 180.0)
        ax.set_theta_direction( -1 )
        ax.set_theta_zero_location( 'W' )

        for iX2 in range( self.nX[1] ):

            self.Data[:,iX2] \
              = ( self.Data[:,iX2])# - self.uK ) / self.uK

        self.vmin = self.Data.min()
        self.vmax = self.Data.max()

        self.Norm = GetNorm( self.UseLogScale, self.Data, \
                             vmin = self.vmin, vmax = self.vmax )

        im = ax.pcolormesh( self.theta, self.r, self.Data, \
                            cmap = self.cmap, \
                            norm = self.Norm, \
                            shading = 'nearest' )

        #ax.set_position( [0.1,-0.45,0.8,2] )
        #cax = fig.add_axes( [0.93,0.1,0.03,0.8] )
        #cbar = fig.colorbar( im, cax = cax )
        cbar = fig.colorbar( im )
        cbar.set_label( self.Field + ' ' + self.DataUnit )
#        cbar.set_label( '(' + self.Field + ' - <' + self.Field + '>)' \
#                        ' / ' + '|<' + self.Field + '>|' )

        ax.set_rmin( 0.0 )
        ax.set_rmax( 150.0 )

        #plt.savefig( 'fig.{:}_{:}.png'.format( ID, Field ), dpi = 300 )
        plt.show()
        plt.close()

if __name__ == '__main__':

    ID = 'GR2D_M2.0_Mdot0.3_Rs150_ColdStart'
    Field = 'AF_P'
    UseLogScale = True

    PlotFields \
      = PlotFieldsAMReX2D( ID, Field, cmap = 'viridis', UseLogScale = UseLogScale )

    PlotFields.GetData()

    PlotFields.ComputeAngleAverage()

    PlotFields.PlotData()

    os.system( 'rm -rf __pycache__' )
