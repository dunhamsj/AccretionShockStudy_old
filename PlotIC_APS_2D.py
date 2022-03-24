#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
from sys import argv

from UtilitiesModule import GetData as GD
from UtilitiesModule import GetNorm

Root = './'

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

        PlotFileBaseName = self.ID + '.plt'

        self.Data, self.DataUnit, self.r, self.theta, \
          self.Time, self.xL, self.xU \
            = GD( self.DataDirectory, PlotFileBaseName, \
                  argv, self.Field )

        self.nX = np.array( [ self.Data.shape[0], self.Data.shape[1], 1 ], \
                            np.int64 )

        self.dX = ( self.xU - self.xL ) / np.float64( self.nX )

        return


    def ComputeAngleAverage( self ):

        self.uK = np.empty( (self.nX[0]), np.float64 )

        theta = self.theta[0]

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

        fig = plt.figure( figsize = (8,6) )
        ax  = fig.add_subplot( 111, polar = True )

#        fig.suptitle( self.ID + '\nTime = ' \
#                      + str( np.int64( self.Time ) ) + ' ms' )

        ax.set_thetamin( 0.0  )
        ax.set_thetamax( 180.0)
        ax.set_theta_direction( 1 )
        ax.set_theta_zero_location( 'N' )

        for iX2 in range( self.nX[1] ):

            self.Data[:,iX2] \
              = ( self.Data[:,iX2] - self.uK ) / self.uK

        self.vmin = self.Data.min()
        self.vmax = self.Data.max()

        self.Norm = GetNorm( self.UseLogScale, self.Data, \
                             vmin = self.vmin, vmax = self.vmax )

        im = ax.pcolormesh( self.theta, self.r, self.Data, \
                            cmap = self.cmap, \
                            norm = self.Norm, \
                            shading = 'nearest' )

        #ax.set_position( [0.1,-0.45,0.8,2] )
        cax = fig.add_axes( [0.75,0.1,0.03,0.8] )
        cbar = fig.colorbar( im, cax = cax )
        cbar.set_label( r'$\Delta p/p\left(r_{c}\right)$', labelpad = -25 )
        mn = self.Data.min()
        mx = self.Data.max()
        cbar.set_ticks( [ mn, mx ] )
        cbar.set_ticklabels( [r'$-10^{-6}$',r'$+10^{-6}$'] )


        ax.set_rmin( 0.0 )
        ax.set_rmax( 185.0 )

        plt.savefig( 'fig.IC_2D_APS.png', dpi = 300, bbox_inches='tight' )
        #plt.show()
        plt.close()

if __name__ == '__main__':

    ID = 'GR2D_M2.0_Mdot0.3_Rs180'
    Field = 'AF_P'
    UseLogScale = True

    PlotFields \
      = PlotFieldsAMReX2D( ID, Field, cmap = 'RdBu', UseLogScale = UseLogScale )

    PlotFields.GetData()

    PlotFields.ComputeAngleAverage()

    PlotFields.PlotData()

    os.system( 'rm -rf __pycache__' )
