#!/usr/bin/env python3

import numpy as np
import os
from multiprocessing import Process, cpu_count

from UtilitiesModule import Overwrite, GetData, GetFileArray

def MakeDataFile \
      ( field, plotFileDirectory, dataFileDirectory, \
        plotFileBaseName, \
        SSi = -1, SSf = -1, nSS = -1, \
        maxLevel = -1, verbose = False, \
        forceChoice = False, OW = True ):

    """
    Generate a directory containing data files where each data file corresponds
    to a specific AMReX plotfile. Mesh data are stored in the header of each
    data file and can be accessed with the ReadHeader function.
    """

    if verbose:
        print( '\n  Running MakeDataFile' )
        print(   '  --------------------' )

    if not dataFileDirectory[-1] == '/' : dataFileDirectory += '/'

    if verbose:
        print( '\n    dataFileDirectory: {:}' \
               .format( dataFileDirectory ) )

    OW = Overwrite( dataFileDirectory, forceChoice, OW )

    if OW:

        os.system( 'rm -rf {:}'.format( dataFileDirectory ) )
        os.system(  'mkdir {:}'.format( dataFileDirectory ) )

        if plotFileDirectory[-1] != '/' : plotFileDirectory += '/'

        if verbose :
            print( '\n    plotFileDirectory: {:}' \
                   .format( plotFileDirectory ) )

        plotFileNameArray = GetFileArray( plotFileDirectory, plotFileBaseName )

        if SSi < 0: SSi = 0
        if SSf < 0: SSf = plotFileNameArray.shape[0] - 1
        if nSS < 0: nSS = plotFileNameArray.shape[0]

        plotFileArray = []
        for i in range( nSS ):
            iSS = SSi + np.int64( ( SSf - SSi ) / ( nSS - 1 ) * i )
            plotFile = str( plotFileNameArray[iSS] )
            if plotFile[-1] == '/' :
                plotFileArray.append( plotFile[0:-1] )
            else:
                plotFileArray.append( plotFile )
        plotFileArray = np.array( plotFileArray )

        timeHeaderBase = '# Time [ms]: '
        X1Base         = '# X1_C [km]: '
        X2Base         = '# X2_C [rad]: '
        X3Base         = '# X3_C [rad]: '
        dX1Base        = '# dX1  [km]: '
        dX2Base        = '# dX2  [rad]: '
        dX3Base        = '# dX3  [rad]: '

        # Just to get dimensionality
        plotFile = plotFileDirectory + plotFileArray[0]
        time, data, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFile, field, verbose )

        nDimsX = 1
        if( nX[1] > 1 ): nDimsX += 1
        if( nX[2] > 1 ): nDimsX += 1

        def loop( iLo, iHi ):

            N = iHi - iLo
            for i in range( iLo, iHi ):

                dataFile = dataFileDirectory + plotFileArray[i] + '.dat'
                plotFile = plotFileDirectory + plotFileArray[i]

                if verbose:
                    print( 'Generating data file: {:} ({:}/{:})'.format \
                             ( dataFile, i+1, N ) )

                time, data, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
                  = GetData( plotFile, field )

                if   nDimsX == 1:
                    dataShape = '{:d}'.format( X1.shape[0] )
                    dshape = (X1.shape[0])
                elif nDimsX == 2:
                    dataShape = '{:d} {:d}'.format( X1.shape[0], X2.shape[0] )
                    dshape = (X1.shape[0],X2.shape[0])
                else:
                    exit( 'MakeDataFile not implemented for nDimsX > 2' )

                data = np.copy( data.reshape( dshape ) )

                # Save multi-D array with np.savetxt. Taken from:
                # https://stackoverflow.com/questions/3685265/
                # how-to-write-a-multidimensional-array-to-a-text-file

                with open( dataFile, 'w' ) as fileOut:

                    fileOut.write( '# {:}\n'             .format( dataFile  ) )
                    fileOut.write( '# Array Shape: {:}\n'.format( dataShape ) )
                    fileOut.write( '# Data Units: {:}\n' .format( dataUnits ) )

                    timeHeader = timeHeaderBase + '{:.16e}\n'.format( time )
                    fileOut.write( timeHeader )

                    fileOut.write( X1Base )
                    for iX1 in range( X1.shape[0] ):
                        fileOut.write( str( X1[iX1] ) + ' ' )
                    fileOut.write( '\n' )

                    fileOut.write( X2Base )
                    for iX2 in range( X2.shape[0] ):
                        fileOut.write( str( X2[iX2] ) + ' ' )
                    fileOut.write( '\n' )

                    fileOut.write( X3Base )
                    for iX3 in range( X3.shape[0] ):
                        fileOut.write( str( X3[iX3] ) + ' ' )
                    fileOut.write( '\n' )

                    fileOut.write( dX1Base )
                    for iX1 in range( dX1.shape[0] ):
                        fileOut.write( str( dX1[iX1] ) + ' ' )
                    fileOut.write( '\n' )

                    fileOut.write( dX2Base )
                    for iX2 in range( dX2.shape[0] ):
                        fileOut.write( str( dX2[iX2] ) + ' ' )
                    fileOut.write( '\n' )

                    fileOut.write( dX3Base )
                    for iX3 in range( dX3.shape[0] ):
                        fileOut.write( str( dX3[iX3] ) + ' ' )
                    fileOut.write( '\n' )

                    np.savetxt( fileOut, data )

                # end with open( DataFileName, 'w' ) as FileOut

            # end for i in range( nSS )

        # end of loop( iLo, iHi )

        # Adapted from:
        # https://www.benmather.info/post/2018-11-24-multiprocessing-in-python/

        nProc = 8#cpu_count()

        print( 'Generating {:} with {:} processes...\n'.format \
             ( dataFileDirectory, nProc ) )

        processes = []

        for i in range( nProc ):
            iLo = np.int64( np.float64( i     ) / np.float64( nProc ) * nSS )
            iHi = np.int64( np.float64( i + 1 ) / np.float64( nProc ) * nSS )
            p = Process( target = loop, args = (iLo,iHi) )
            p.start()
            processes.append( p )

        [ p.join() for p in processes ]

    # end if OW

    return # end of MakeDataFile

def ReadHeader( DataFile ):

    f = open( DataFile )

    dum = f.readline()

    s = f.readline(); ind = s.find( ':' )+1
    dataShape = np.array( list( map( np.int64, s[ind:].split() ) ), np.int64 )

    s = f.readline(); ind = s.find( ':' )+1
    dataUnits = s[ind:]

    s = f.readline(); ind = s.find( ':' )+1
    time = np.float64( s[ind:] )

    s = f.readline(); ind = s.find( ':' )+1
    X1_C = np.array( list( map( np.float64, s[ind:].split() ) ), np.float64 )

    s = f.readline(); ind = s.find( ':' )+1
    X2_C = np.array( list( map( np.float64, s[ind:].split() ) ), np.float64 )

    s = f.readline(); ind = s.find( ':' )+1
    X3_C = np.array( list( map( np.float64, s[ind:].split() ) ), np.float64 )

    s = f.readline(); ind = s.find( ':' )+1
    dX1 = np.array( list( map( np.float64, s[ind:].split() ) ), np.float64 )

    s = f.readline(); ind = s.find( ':' )+1
    dX2 = np.array( list( map( np.float64, s[ind:].split() ) ), np.float64 )

    s = f.readline(); ind = s.find( ':' )+1
    dX3 = np.array( list( map( np.float64, s[ind:].split() ) ), np.float64 )

    f.close()

    return dataShape, dataUnits, time, X1_C, X2_C, X3_C, dX1, dX2, dX3
