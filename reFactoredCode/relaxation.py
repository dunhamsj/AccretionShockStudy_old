#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from multiprocessing import Process, cpu_count, Manager

from UtilitiesModule import Overwrite, GetFileArray, GetData

def getData( plotfileDirectory, ID, field, nX, forceChoice, OW ):

    dataFileName = '.{:}_Relaxation_{:}_nX{:}.dat' \
                   .format( ID, field, str( nX ).zfill(4) )

    OW = Overwrite( dataFileName, ForceChoice = forceChoice, OW = OW )

    if OW:

        plotfileBaseName = ID + '.plt'

        plotfileArray = GetFileArray( plotfileDirectory, plotfileBaseName )
        plotfileArray = np.copy( plotfileArray[:-1] )

        nSS = plotfileArray.shape[0]

        ind = np.linspace( 0, nSS-1, nSS, dtype = np.int64 )

        print()
        print( '  Generating derivative data' )
        print( '  --------------------------' )

        Data     = np.zeros( (nSS,nX), np.float64 )
        Gradient = np.zeros( (nSS-1) , np.float64 )
        Time     = np.zeros( (nSS)   , np.float64 )

        def loop( iLo, iHi, nodalData, time, iProc, \
                  return_data, return_time ):

            for i in range( iLo, iHi ):

                N = iHi - iLo
                j = i - iLo

                if j % 10 == 0:
                    print( '    File {:d}/{:d}'.format( j+1, N ) )

                data, dataUnits, time[j] \
                  = GetData( plotfileDirectory, plotfileBaseName, field, \
                             'spherical', True, argv = ['a',plotfileArray[i]], \
                             ReturnTime = True, ReturnMesh = False, \
                             Verbose = False )

                nodalData[j] = np.copy( data[:,0,0] )

            return_data[iProc] = nodalData
            return_time[iProc] = time

        # Adapted from:
        # https://www.benmather.info/post/
        # 2018-11-24-multiprocessing-in-python/

        nProcs = max( 10, cpu_count() // 2 )

        manager     = Manager()
        return_data = manager.dict()
        return_time = manager.dict()

        if nProcs > 1:

            processes = []

            for iProc in range( nProcs ):

                iLo \
                  = np.int64( np.float64( iProc     ) \
                      / np.float64( nProcs ) * nSS )
                iHi \
                  = np.int64( np.float64( iProc + 1 ) \
                      / np.float64( nProcs ) * nSS )

                d = np.empty( (iHi-iLo,nX), np.float64 )
                t = np.empty( (iHi-iLo)   , np.float64 )

                p = Process \
                      ( target = loop, \
                        args = (iLo,iHi,d,t,iProc,return_data,return_time) )
                p.start()
                processes.append( p )

            # MPI BARRIER
            [ p.join() for p in processes ]

            for iProc in range( nProcs ):

                iLo \
                  = np.int64( np.float64( iProc     ) \
                      / np.float64( nProcs ) * nSS )
                iHi \
                  = np.int64( np.float64( iProc + 1 ) \
                      / np.float64( nProcs ) * nSS )

                Data[iLo:iHi] = return_data[iProc]
                Time[iLo:iHi] = return_time[iProc]

        else:

            loop( 0, nSS-1, Data, Time, 0, return_data, return_time )

        for i in range( Gradient.shape[0] ):

            Num = ( Data[i] - Data[i-1] ) / ( 0.5 * ( Data[i] + Data[i-1] ) )
            Den = Time[i+1] - Time[i-1]

            Gradient[i-1] = ( Num / Den ).max()

        np.savetxt( dataFileName, np.vstack( (Time[:-1],Gradient) ) )

        del plotfileBaseName, plotfileArray, \
            Data, Gradient, Time

    # END if OW

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

    nX = 460

    UseLogScale = True

    #ID = 'NR1D_M1.4_Rpns040_Rs1.80e2_ST1.00e-6'
    ID = 'NR1D_M1.4_Rpns040_Rs1.80e2_ST1.00e-6'

    SaveFileAs = 'fig.Relaxation_{:}.png'.format( ID )

    Root = '/lump/data/accretionShockStudy/newData/'

    plotfileDirectory = Root + ID + '/'
    ID = 'NR1D_M1.4_Rpns040_Rs180_Mdot0.3'

    D = 'PF_D'
    Time, Data = getData( plotfileDirectory, ID, D, nX, False, True )

    ind = np.where( Time >= 0.0 )[0]

    Time = np.copy( Time[ind] )

    fig, ax = plt.subplots( 1, 1 )

    fig.suptitle( ID )

    PlotRelaxationVsTime( ax, Time, Data, D, ID, UseLogScale )

    ax.grid()

#    plt.show()

    plt.savefig( SaveFileAs, dpi = 300 )

    plt.close()

    import os
    os.system( 'rm -rf __pycache__ ' )
