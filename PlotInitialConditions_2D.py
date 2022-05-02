#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
from sys import argv

from UtilitiesModule import GetData as GD
from UtilitiesModule import GetNorm

Root = '/lump/data/AccretionShockStudy/'

class PlotFieldsAMReX2D:

    def __init__( self, ID, Field, UseLogScale = True, \
                  suffix = '' ):

        self.ID = ID
        self.Field = Field

        self.DataDirectory = Root + ID + '{:}/'.format( suffix )
        self.UseLogScale = UseLogScale

        return

    def GetData( self ):

        PlotFileBaseName = self.ID + '.plt'

        self.Data, self.DataUnit, self.r, self.theta, \
          self.Time, self.xL, self.xU \
            = GD( self.DataDirectory, PlotFileBaseName, \
                  self.Field, argv = [ 'a', '0' ] )

        self.r     = np.copy( self.r    [:,0] )
        self.theta = np.copy( self.theta[0,:] )

        self.nX = np.array( [ self.r.shape[0], self.theta.shape[0], 1 ], \
                            np.int64 )

        self.dX = ( self.xU - self.xL ) / np.float64( self.nX )

        self.vmin = self.Data.min()
        self.vmax = self.Data.max()

        return


    def ComputeAngleAverage( self ):

        self.uK = np.empty( (self.nX[0]), np.float64 )

        for iX1 in range( self.nX[0] ):

            self.uK[iX1] \
              = 0.5 * np.sum( self.Data[iX1] * np.sin( self.theta ) ) \
                  * self.dX[1]

        return

    def PlotData( self ):

        fig = plt.figure( figsize = (8,6) )
        ax  = fig.add_subplot( 111 )

        fig.suptitle( self.ID + '\nTime = ' \
                      + str( np.int64( self.Time ) ) + ' ms' )

#        for iX2 in range( self.nX[1] ):
#
#            self.Data[:,iX2] \
#              = ( self.Data[:,iX2] - self.uK ) / self.uK

        self.Norm = GetNorm( self.UseLogScale, self.Data, \
                             vmin = self.vmin, vmax = self.vmax )

        iX2 = 0
        im = ax.plot( self.r, self.Data[:,iX2] )

        ax.set_xlim( self.xL[0], self.xU[0] )

        #plt.savefig( 'fig.IC_2D.png', dpi = 300, bbox_inches='tight' )
        plt.show()
        plt.close()

if __name__ == '__main__':

    IDs = np.array( [ 'GR2D_M2.0_Mdot0.1_Rs150', \
                      'GR2D_M2.0_Mdot0.3_Rs150', \
                      'GR2D_M2.0_Mdot0.5_Rs150' ] )

    Fields      = np.array( [ 'AF_P' ] )
    UseLogScale = True

    dum = PlotFieldsAMReX2D( IDs[0], 'PF_D', UseLogScale = False )
    dum.GetData()
    Time = dum.Time
    xL   = dum.xL
    xU   = dum.xU
    r    = dum.r

    fig = plt.figure( figsize = (9,6) )
    ax  = fig.add_subplot( 111 )

    fig.suptitle( IDs[0][0:10]+IDs[0][-5:]+'\nTime = ' \
                  + str( np.int64( Time ) ) + ' ms' )

    c  = [ 'r', 'b', 'm' ] # one for each ID
    ls = [ '-', '--' ] # one for each field

    for id in range( IDs.shape[0] ):

        for field in range( Fields.shape[0] ):

            d \
              = PlotFieldsAMReX2D( IDs[id], Fields[field], \
                                   UseLogScale = UseLogScale )
            d.GetData()

#            Data.ComputeAngleAverage()

            iX2 = 0
            ax.plot( r, d.Data[:,iX2], label = IDs[id][10:17] + '_{:}'.format \
                     ( Fields[field] ), c = c[id], ls = ls[field] )

    if UseLogScale: ax.set_yscale( 'log' )
    ax.set_ylabel( 'erg/cm^3', fontsize = 15 )
    ax.set_xlabel( 'r [km]', fontsize = 15 )
    ax.tick_params( axis='both', which='major', labelsize=13 )
    ax.grid()

    ax.legend(prop={'size':13})
    ax.set_xlim( xL[0], xU[0] )

    plt.savefig( 'fig.IC_2D_P.png', dpi = 300, bbox_inches='tight' )
    #plt.show()
    plt.close()

    os.system( 'rm -rf __pycache__' )
