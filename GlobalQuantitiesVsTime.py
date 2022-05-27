#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
from sys import argv
plt.style.use( './Publication.sty' )

from UtilitiesModule import GetData as GD, GetFileArray, Overwrite

class GlobalQuantities:

    def __init__( self, Root, ID ):

        self.Root = Root
        self.ID = ID

        self.DataDirectory = self.Root + self.ID + '/'
        self.PlotFileBaseName = self.ID + '.plt'

        return

    def IntegrateField( self, q, SqrtGm, dX1, dX2 ):

        nX1 = dX1.shape[0]
        nX2 = dX2.shape[0]

        Q = 0.0
        for iX1 in range( nX1 ):
            for iX2 in range( nX2 ):
                Q += 2.0 * np.pi * dX1[iX1] * dX2[iX2] \
                       * q[iX1,iX2] * SqrtGm[iX1,iX2]

        return Q

    def GetSnapshot( self, Field, iSS ):

        Data, DataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xU, nX, Time \
          = GD( self.DataDirectory, self.PlotFileBaseName, Field, \
                'spherical', True, argv = [ 'a', iSS ], \
                ReturnTime = True, ReturnMesh = True, Verbose = False )

        SqrtGm, DataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xU, nX, Time \
          = GD( self.DataDirectory, self.PlotFileBaseName, 'GF_SqrtGm', \
                'spherical', True, argv = [ 'a', iSS ], \
                ReturnTime = True, ReturnMesh = True, Verbose = False )

        SqrtGm *= 1.0e5**2

        return Data, SqrtGm, dX1, dX2, Time

    def GetData( self, Field ):

        FileName = self.ID + '_' + Field + '.dat'

        OW = Overwrite( FileName, ForceChoice = True, OW = False )

        if OW:

            FileArray \
              = GetFileArray( self.DataDirectory, self.PlotFileBaseName )#[::10]

            Q    = np.empty( FileArray.shape[0], np.float64 )
            Time = np.empty( FileArray.shape[0], np.float64 )

            for iSS in range( FileArray.shape[0] ):

                print( '{:d}/{:d}'.format( iSS+1, FileArray.shape[0] ) )

                Data, SqrtGm, dX1, dX2, Time[iSS] \
                  = self.GetSnapshot( Field, FileArray[iSS] )

                Q[iSS] = self.IntegrateField( Data, SqrtGm, dX1, dX2 )

            np.savetxt( FileName, np.vstack( ( Time, Q ) ) )

        d = np.loadtxt( FileName )

        Time = d[0]
        Q    = d[1]

        return Time, Q

if __name__ == '__main__':

    Root = '/lump/data/AccretionShockStudy/'
    IDs = np.array( [ 'GR2D_M2.0_Mdot0.3_Rs150', \
                      'NR2D_M2.0_Mdot0.3_Rs150' ], str )

    Fields = np.array( [ 'PF_D', 'TurbulentEnergyDensity' ], str )

    GR = GlobalQuantities( Root, IDs[0] )
    NR = GlobalQuantities( Root, IDs[1] )

    Time, GR_Mass = GR.GetData( Fields[0] )
    Time, GR_TE   = GR.GetData( Fields[1] )
    Time, NR_Mass = NR.GetData( Fields[0] )
    Time, NR_TE   = NR.GetData( Fields[1] )

    fig, axs = plt.subplots( Fields.shape[0], 1 )
    axs[-1].set_xlabel( 'Time [ms]' )

    c = [ 'r', 'b' ]
    #axs[0].plot( Time, ( GR_Mass - GR_Mass[0] ) / GR_Mass[0], \
    #             c = c[0], label = 'GR' )
    #axs[0].plot( Time, ( NR_Mass - NR_Mass[0] ) / NR_Mass[0], \
    #             c = c[1], label = 'NR' )
    #axs[1].plot( Time, ( GR_TE   - GR_TE[0]   ) / GR_TE[0]  , \
    #             c = c[0] )
    #axs[1].plot( Time, ( NR_TE   - NR_TE[0]   ) / NR_TE[0]  , \
    #             c = c[1] )
    axs[0].plot( Time, GR_Mass, \
                 c = c[0], label = 'GR' )
    axs[0].plot( Time, NR_Mass, \
                 c = c[1], label = 'NR' )
    axs[1].plot( Time, GR_TE  , \
                 c = c[0] )
    axs[1].plot( Time, NR_TE , \
                 c = c[1] )

    axs[0].legend()
    axs[0].set_ylabel( 'M(t)' )
    axs[1].set_ylabel( 'TE(t)' )
    axs[0].set_yscale( 'log' )
    axs[1].set_yscale( 'log' )
    axs[1].set_ylim( 1.0e32 )

    plt.show()

    os.system( 'rm -rf __pycache__' )
