#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from matplotlib import animation
from multiprocessing import Process, cpu_count, Manager

class Relaxation:

    #Root = '/home/dunhamsj/AccretionShockData/'
    Root = '/lump/data/accretionShockStudy/'

    def __init__( self, nX, SSi, SSf, nSS = -1, \
                  ForceChoice = False, Overwrite = False ):

        self.nX  = nX
        self.SSi = SSi
        self.SSf = SSf
        self.nSS = nSS

        if self.nSS < 0: self.nSS = SSf - SSi

        print( '' )
        print( '  Creating instance of Relaxation...' )
        print( '  ----------------------------------' )
        print( '{:>8} {:}'.format( 'nX:' , self.nX ) )
        print( '{:>8} {:}'.format( 'SSi:', self.SSi ) )
        print( '{:>8} {:}'.format( 'SSf:', self.SSf ) )
        print( '{:>8} {:}'.format( 'nSS:', self.nSS ) )
        print( '' )

        self.ind = np.linspace( self.SSi, self.SSf, self.nSS, dtype = np.int64 )

        self.ForceChoice = ForceChoice
        self.OW = Overwrite

        return


    def SetDataDirectory( self, ID, suffix = '' ):

        self.DataDirectory = self.Root + ID + '{:}/'.format( suffix )

        print( '{:>18} {:}'.format( 'ID:', ID ) )
        print( '{:>18} {:}'.format( 'DataDirectory:', self.DataDirectory ) )
        print( '' )

        return


    def GetData( self, Field, ID ):

        import UtilitiesModule as UM

        self.SetDataDirectory( ID )

        print( '{:>10} {:}'.format( 'Field:', Field ) )

        FileName = '.{:}_Relaxation_{:}.dat'.format( ID, Field )

        OW = UM.Overwrite( FileName, self.ForceChoice, self.OW )

        if OW:

            print( '    Generating File: {:}'.format( FileName ) )

            PlotFileBaseName = ID + '.plt'

            FileArray = UM.GetFileArray( self.DataDirectory, PlotFileBaseName )
            File \
              = UM.ChoosePlotFile( self.DataDirectory, PlotFileBaseName, argv, \
                                   Verbose = False )

            def loop( iLo, iHi, NodalData, Time, iProc, \
                      return_data, return_time ):

                for i in range( iLo, iHi ):

                    j = i - iLo

                    if i % 10 == 0:
                        print( 'File {:d}/{:d}'.format( i, self.nSS ) )

                    NodalData[j], DataUnit, r, theta, phi, dr, dtheta, dphi, \
                      xL, xU, nX, Time[j] \
                        = UM.GetData( self.DataDirectory, PlotFileBaseName, \
                                      Field, 'spherical', True, \
                                      argv = [ 'a', FileArray[self.ind[i]] ], \
                                      ReturnTime = True, ReturnMesh = True, \
                                      Verbose = False )

                return_data[iProc] = NodalData
                return_time[iProc] = Time

            # Adapted from:
            # https://www.benmather.info/post/
            # 2018-11-24-multiprocessing-in-python/

            nProcs = cpu_count()

            processes = []

            nSS = self.ind.shape[0]

            manager = Manager()
            return_data = manager.dict()
            return_time = manager.dict()

            for i in range( nProcs ):

                iLo \
                  = np.int64( np.float64( i     ) / np.float64( nProcs ) * nSS )
                iHi \
                  = np.int64( np.float64( i + 1 ) / np.float64( nProcs ) * nSS )

                nd = np.empty( (iHi-iLo,self.nX), np.float64 )
                t  = np.empty( (iHi-iLo)        , np.float64 )

                p = Process \
                      ( target = loop, \
                        args = (iLo,iHi,nd,t,i,return_data,return_time) )
                p.start()
                processes.append( p )

            # MPI BARRIER
            [ p.join() for p in processes ]

            NodalData = np.zeros( (self.nSS,self.nX), np.float64 )
            DiffData  = np.zeros( (self.nSS-1)      , np.float64 )
            Time      = np.zeros( (self.nSS)        , np.float64 )

            for i in range( nProcs ):
                iLo \
                  = np.int64( np.float64( i     ) / np.float64( nProcs ) * nSS )
                iHi \
                  = np.int64( np.float64( i + 1 ) / np.float64( nProcs ) * nSS )

                NodalData[iLo:iHi] = return_data[i]
                Time     [iLo:iHi] = return_time[i]

            Den = np.empty( (self.nX), np.float64 )

            for i in range( DiffData.shape[0] ):

                Num = np.abs( NodalData[i+1] - NodalData[i] )
                for j in range( Num.shape[0] ):

                    Den[j] \
                      = max( 1.0e-17, \
                             0.5 * np.abs( NodalData[i+1,j] + NodalData[i,j] ) )

                DiffData[i] = ( Num / Den ).max()

            np.savetxt( FileName, np.vstack( (Time[:-1],DiffData) ) )

            del PlotFileBaseName, FileArray, File, \
                NodalData, DiffData, Time, Num, Den

        Time, Data = np.loadtxt( FileName )

        return Time, Data


    def PlotRelaxationVsTime \
          ( self, ax, Time, Data, Field, ID, UseLogScale, label = '' ):

        ax.plot( Time, np.abs( Data ), '.', \
                 markersize = 2.0, markevery = 1, label = label )

        ax.set_xlabel( 'Time [ms]' )
        ax.set_ylabel( 'max( |d{:s}/dt / {:s}| )'.format( Field, Field ) )

        if UseLogScale: ax.set_yscale( 'log' )

        return

if __name__ == '__main__':

    nX  = 640
    SSi = 0
    SSf = 1999
    nSS = 2000

    Relax \
      = Relaxation( nX, SSi, SSf, nSS, \
                    ForceChoice = False, Overwrite = True )

    UseLogScale = False

    ID = 'GR1D_M2.0_Mdot0.3_Rs150_Gm1.13'

    SaveFileAs = 'fig.Relaxation_{:}.png'.format( ID )

    D = 'PF_D'
    Time_D, Data_D = Relax.GetData( D, ID )

    V = 'PF_V1'
    Time_V, Data_V = Relax.GetData( V, ID )

    P = 'AF_P'
    Time_P, Data_P = Relax.GetData( P, ID )

    fig, axs = plt.subplots( 3, 1 )

#    fig.suptitle( 'Gaussian perturbation below shock\n{:}'.format( ID ) )

    Relax.PlotRelaxationVsTime( axs[0], Time_D, Data_D, D, ID, UseLogScale )
    Relax.PlotRelaxationVsTime( axs[1], Time_V, Data_V, V, ID, UseLogScale )
    Relax.PlotRelaxationVsTime( axs[2], Time_P, Data_P, P, ID, UseLogScale )

#    axs[0].set_ylim( 1.0e-6, 5.0 )
#    axs[1].set_ylim( 1.0e-6, 5.0 )
#    axs[2].set_ylim( 1.0e-6, 5.0 )
    axs[0].grid()
    axs[1].grid()
    axs[2].grid()

#    plt.show()

    plt.savefig( SaveFileAs, dpi = 300 )

    plt.close()

    import os
    os.system( 'rm -rf __pycache__ ' )
