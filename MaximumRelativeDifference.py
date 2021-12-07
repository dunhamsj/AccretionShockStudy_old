#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from matplotlib import animation

class MaxRelDiff:

    Root = '/Users/dunhamsj/Research/Data/AccretionShockParameterStudy/'

    def __init__( self, nX, SSi, SSf, nSS = -1, \
                  ForceOverwrite = False, Overwrite = False ):

        self.nX  = nX
        self.SSi = SSi
        self.SSf = SSf
        self.nSS = nSS

        if self.nSS < 0: self.nSS = SSf - SSi

        print( '' )
        print( '  Creating instance of MaxRelDiff...' )
        print( '  ----------------------------------' )
        print( '    nX:  {:}'.format( self.nX ) )
        print( '    SSi: {:}'.format( self.SSi ) )
        print( '    SSf: {:}'.format( self.SSf ) )
        print( '    nSS: {:}'.format( self.nSS ) )
        print( '' )

        self.ind = np.linspace( self.SSi, self.SSf, self.nSS, dtype = np.int64 )

        self.ForceOverwrite = ForceOverwrite
        self.Overwrite = Overwrite

        return


    def SetDataDirectory( self, ID, suffix = '' ):

        self.DataDirectory = self.Root + ID + '{:}/'.format( suffix )

        print( '    ID:            {:}'.format( ID ) )
        print( '    DataDirectory: {:}'.format( self.DataDirectory ) )
        print( '' )

        return


    def GetData( self, Field, ID ):

        import UtilitiesModule as UM

        self.SetDataDirectory( ID )

        FileName = '.{:}_MaxRelDiff_{:}.dat'.format( ID, Field )

        OW = UM.OverwriteFile( FileName, self.ForceOverwrite, self.Overwrite )

        if OW:

            PlotFileBaseName = ID + '.plt'

            FileArray = UM.GetFileArray( self.DataDirectory, PlotFileBaseName )
            File \
              = UM.ChoosePlotFile( self.DataDirectory, PlotFileBaseName, argv, \
                                   Verbose = False )

            NodalData = np.empty( (self.nSS,self.nX), np.float64 )
            DiffData  = np.empty( (self.nSS-1), np.float64 )
            Time      = np.empty( (self.nSS), np.float64 )

            for i in range( self.ind.shape[0] ):

                if i % 10 == 0:
                    print( 'File {:d}/{:d}'.format( i, self.nSS ) )

                NodalData[i], DataUnit, r, theta, Time[i], xL, xU \
                  = UM.GetData( self.DataDirectory, PlotFileBaseName, \
                                [ 'a', FileArray[self.ind[i]] ], Field, \
                                Verbose = False )

            Den = np.empty( (self.nX), np.float64 )

            for i in range( DiffData.shape[0] ):

                Num = np.abs( NodalData[i+1] - NodalData[i] )
                for j in range( Num.shape[0] ):

                    Den[j] \
                      = max( 1.0e-17, \
                             0.5 * np.abs( NodalData[i+1,j] + NodalData[i,j] ) )

                DiffData[i] = ( Num / Den ).max()

            np.savetxt( FileName, np.vstack( (Time[:-1],DiffData) ) )

        Time, Data = np.loadtxt( FileName )

        return Time, Data


    def PlotMaxRelDiffVsTime \
          ( self, ax, Time, Data, Field, ID, UseLogScale, label = '' ):

        ax.plot( Time, np.abs( Data ), '.', \
                 markersize = 2.0, markevery = 1, label = label )

        ax.set_xlabel( 'Time [ms]' )
        ax.set_ylabel( 'MRD ({:s})'.format( Field ) )

        if UseLogScale: ax.set_yscale( 'log' )

        return

if __name__ == '__main__':

    nX  = 640
    SSi = 0
    SSf = 1999
    nSS = 2000

    MRD = MaxRelDiff( nX, SSi, SSf, nSS )

    UseLogScale = True

    Field = 'PF_D'

    ID = 'NR1D_M1.4_Mdot0.3_Rs180_PA1.00e-05_nX640'

    SaveFileAs = 'fig.MaxRelDiff_{:}.png'.format( ID )

    D = 'PF_D'
    Time_D, Data_D = MRD.GetData( D, ID )

    V = 'PF_V1'
    Time_V, Data_V = MRD.GetData( V, ID )

    P = 'AF_P'
    Time_P, Data_P = MRD.GetData( P, ID )

    fig, axs = plt.subplots( 3, 1 )

    fig.suptitle( 'Gaussian perturbation below shock\n{:}'.format( ID ) )

    MRD.PlotMaxRelDiffVsTime( axs[0], Time_D, Data_D, D, ID, UseLogScale )
    MRD.PlotMaxRelDiffVsTime( axs[1], Time_V, Data_V, V, ID, UseLogScale )
    MRD.PlotMaxRelDiffVsTime( axs[2], Time_P, Data_P, P, ID, UseLogScale )

    axs[0].set_ylim( 1.0e-6, 5.0 )
    axs[1].set_ylim( 1.0e-6, 5.0 )
    axs[2].set_ylim( 1.0e-6, 5.0 )

#    plt.show()

    plt.savefig( SaveFileAs, dpi = 300 )

    plt.close()

    import os
    os.system( 'rm -rf __pycache__ ' )
