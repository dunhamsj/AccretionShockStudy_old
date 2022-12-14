#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from multiprocessing import Process, cpu_count, Manager

from UtilitiesModule import Overwrite, GetFileArray, GetData

def getData( plotFileDirectory, ID, field, nX, forceChoice, OW ):

    dataFileName = '.{:}_Relaxation_{:}_nX{:}.dat' \
                   .format( ID, field, str( nX ).zfill(4) )

    OW = Overwrite( dataFileName, ForceChoice = forceChoice, OW = OW )

    if OW:

        plotFileBaseName = ID + '.plt'

        plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )

        nSS = plotFileArray.shape[0]

        ind = np.linspace( 0, nSS-1, nSS, dtype = np.int64 )

        def loop( iLo, iHi, nodalData, Time, iProc, \
                  return_data, return_time ):

            for i in range( iLo, iHi ):

                N = iHi - iLo
                j = i - iLo

                if j % 10 == 0:
                    print( 'File {:d}/{:d}'.format( j, N ) )

                plotFile = plotFileDirectory + plotFileArray[i]

                Time[j], data, dataUnits, r, theta, phi, \
                dr, dtheta, dphi, nX \
                  = GetData( plotFile, field )

                nodalData[j] = np.copy( data[:,0,0] )

#            return_data[iProc] = nodalData
#            return_time[iProc] = Time

#        # Adapted from:
#        # https://www.benmather.info/post/
#        # 2018-11-24-multiprocessing-in-python/
#
#        nProcs = min( 8, cpu_count() )
#
#        processes = []
#
#        manager     = Manager()
#        return_data = manager.dict()
#        return_time = manager.dict()
#
#        for iProc in range( nProcs ):
#
#            iLo \
#              = np.int64( np.float64( iProc     ) / np.float64( nProcs ) * nSS )
#            iHi \
#              = np.int64( np.float64( iProc + 1 ) / np.float64( nProcs ) * nSS )
#
#            nd = np.empty( (iHi-iLo,nX), np.float64 )
#            t  = np.empty( (iHi-iLo)   , np.float64 )
#
#            p = Process \
#                  ( target = loop, \
#                    args = (iLo,iHi,nd,t,iProc,return_data,return_time) )
#            p.start()
#            processes.append( p )
#
#        # MPI BARRIER
#        [ p.join() for p in processes ]
#
#        nodalData = np.zeros( (nSS,nX), np.float64 )
#        DiffData  = np.zeros( (nSS-1)      , np.float64 )
#        Time      = np.zeros( (nSS)        , np.float64 )
#
#        for iProc in range( nProcs ):
#
#            iLo \
#              = np.int64( np.float64( iProc     ) / np.float64( nProcs ) * nSS )
#            iHi \
#              = np.int64( np.float64( iProc + 1 ) / np.float64( nProcs ) * nSS )
#
#            nodalData[iLo:iHi] = return_data[iProc]
#            Time     [iLo:iHi] = return_time[iProc]

        nodalData = np.zeros( (nSS,nX), np.float64 )
        DiffData  = np.zeros( (nSS-1)      , np.float64 )
        Time      = np.zeros( (nSS)        , np.float64 )

        loop( 0, nSS-1, nodalData, Time, 0, 0, 0 )
        Den = np.empty( (nX), np.float64 )

        for i in range( DiffData.shape[0] ):

            Num = np.abs( nodalData[i+1] - nodalData[i] )
            for j in range( Num.shape[0] ):

                Den[j] \
                  = max( 1.0e-17, \
                         0.5 * np.abs( nodalData[i+1,j] + nodalData[i,j] ) )

            DiffData[i] = ( Num / Den ).max()

        np.savetxt( dataFileName, np.vstack( (Time[:-1],DiffData) ) )

        del plotFileBaseName, plotFileArray, \
            nodalData, DiffData, Time, Den

    Time, Data = np.loadtxt( dataFileName )

    return Time, Data


def PlotRelaxationVsTime \
      ( ax, Time, Data, field, ID, UseLogScale, label = '' ):

    ax.plot( Time, np.abs( Data ), '.', \
             markersize = 2.0, markevery = 1, label = label )

    ax.set_xlabel( 'Time [ms]' )

    if field == 'PF_D':
        ylabel = r'max($\dot{\rho}/\rho$)'
    elif field == 'PF_V1':
        ylabel = r'max($\dot{v}/v$)'
    elif field == 'AF_P':
        ylabel = r'max($\dot{p}/p$)'
    else:
        ylabel = ''

    ax.set_ylabel( ylabel )

    if UseLogScale: ax.set_yscale( 'log' )

    return

if __name__ == '__main__':

#    nX = 128
#    nX = 256
#    nX = 384
    nX = 512

    UseLogScale = True

    ID = 'GR1D_M2.8_Mdot0.3_Rs9.00e1_RPNS2.00e1'

    SaveFileAs = 'fig.Relaxation_{:}.png'.format( ID )

#    Root = '/home/dunhamsj/AccretionShockData/'
    Root = '/lump/data/accretionShockStudy/'

    plotFileDirectory = Root + ID + '/'

    D = 'PF_D'
    Time_D, Data_D = getData( plotFileDirectory, ID, D, nX, False, True )

    V = 'PF_V1'
    Time_V, Data_V = getData( plotFileDirectory, ID, V, nX, False, True )

    P = 'AF_P'
    Time_P, Data_P = getData( plotFileDirectory, ID, P, nX, False, True )

    ind = np.where( Time_D >= 0.0 )[0]

    Time_D = np.copy( Time_D[ind] )
    Time_V = np.copy( Time_V[ind] )
    Time_P = np.copy( Time_P[ind] )
    Data_D = np.copy( Data_D[ind] )
    Data_V = np.copy( Data_V[ind] )
    Data_P = np.copy( Data_P[ind] )

    fig, axs = plt.subplots( 3, 1 )

    fig.suptitle( ID )

    PlotRelaxationVsTime( axs[0], Time_D, Data_D, D, ID, UseLogScale )
    PlotRelaxationVsTime( axs[1], Time_V, Data_V, V, ID, UseLogScale )
    PlotRelaxationVsTime( axs[2], Time_P, Data_P, P, ID, UseLogScale )

    for i in range( axs.shape[0] ):
        axs[i].grid()

#    plt.show()

    plt.savefig( SaveFileAs, dpi = 300 )

    plt.close()

    import os
    os.system( 'rm -rf __pycache__ ' )
