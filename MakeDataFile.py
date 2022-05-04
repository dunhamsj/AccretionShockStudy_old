#!/usr/bin/env python3

import yt
import numpy as np
from os import listdir
from os.path import isfile
from sys import exit

from UtilitiesModule import GetData, Overwrite

def MakeDataFile( Field, DataDirectory, DataFileName, \
                  PlotFileBaseName, UsePhysicalUnits = True, \
                  WriteExtras = False, SSi = -1, SSf = -1, nSS = -1 ):

    print( '\nRunning MakeDataFile...\n' )

    print( 'DataDirectory: {:}\n'.format( DataDirectory ) )

    if( not DataFileName[0] == '.' ): DataFileName = '.' + DataFileName

    yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

    c = 2.99792458e10
    if( not UsePhysicalUnits ): c = 1.0

    # Get array of all plot-files

    # Append "/" to DataDirectory, if not present
    if( not DataDirectory[-1] == '/' ): DataDirectory += '/'

    FileArray \
      = np.sort( np.array( [ file for file in listdir( DataDirectory ) ] ) )

    FileList = []

    for iFile in range( FileArray.shape[0] ):

        sFile = FileArray[iFile]

        if( sFile[0:len(PlotFileBaseName)+1] == PlotFileBaseName + '_' \
              and sFile[len(PlotFileBaseName)+1].isdigit() ):
            FileList.append( sFile )

    FileArray = np.array( FileList )

    # Get some general info about the computational domain

    ds       = yt.load( DataDirectory + FileArray[0] )
    MaxLevel = ds.index.max_level
    nX       = ds.domain_dimensions
    xL       = ds.domain_left_edge.to_ndarray()
    xU       = ds.domain_right_edge.to_ndarray()

    nDimsX = 1
    if( nX[1] > 1 ): nDimsX += 1
    if( nX[2] > 1 ): nDimsX += 1

    if( nDimsX == 3 ):

      print( 'FATAL ERROR: MakeDataFile not implemented for nDimsX = 3' )
      exit( 'Exiting...' )

    if( WriteExtras ):

        # This is if you're running on a computing cluster and don't want
        # to copy the files to your local machine

        with open( 'FileArray.txt', 'w' ) as f:
            for i in FileArray:
                f.write( i )
                f.write( '\n' )
        with open( 'Numbers.txt', 'w' ) as f:
            for i in nX:
                f.write( str(i) )
                f.write( ' ' )
            f.write( '\n' )
            for i in xL.to_ndarray():
                f.write( str(i) )
                f.write( ' ' )
            f.write( '\n' )
            for i in xU.to_ndarray():
                f.write( str(i) )
                f.write( ' ' )
            f.write( '\n' )

        exit()

    if SSi < 0: SSi = 0
    if SSf < 0: SSf = FileArray.shape[0]-1
    if nSS < 0: nSS = SSf - SSi + 1

    GenerateFile = Overwrite( DataFileName )

#    FileList = []
#    if GenerateFile:
#        for i in range( SSi, SSf ):
#            ds = yt.load( '{:}'.format( DataDirectory + FileArray[i] ) )
#    else:
#        for i in range( SSi, SSf ):
#            FileList.append( 'a' )
#        FileArray = np.array( FileList )

    if( GenerateFile ):

        # Put all time-slices into one array to use for movie making

        if  ( nDimsX == 1 ): DataShape = (nSS,nX[0])
        elif( nDimsX == 2 ): DataShape = (nSS,nX[0],nX[1])

        Data = np.empty( DataShape, np.float64 )
        Time = np.empty( nSS, np.float64 )

        print( 'Generating data file: {:}...'.format( DataFileName ) )

        for i in range( nSS ):

            iSS = SSi + np.int64( ( SSf - SSi ) / nSS ) * i
            if i % 10 == 0:
                print( '{:}/{:}'.format( i, nSS ) )

            ds = yt.load( '{:}'.format( DataDirectory + FileArray[iSS] ) )

            Data[i], DataUnit, r, theta, tt, xx, xl \
              = GetData( DataDirectory, PlotFileBaseName, Field, \
                         [ 'a', FileArray[iSS] ] )

            Time[i] = ds.current_time

        # Save multi-D array with np.savetxt. Taken from:
        # https://stackoverflow.com/questions/3685265/
        # how-to-write-a-multidimensional-array-to-a-text-file

        with open( DataFileName, 'w' ) as FileOut:

            FileOut.write( '# Array shape: {:}\n'.format( DataShape ) )
            FileOut.write( '# Units: {:}\n'.format( DataUnit ) )

        with open( DataFileName, 'a' ) as FileOut:

            FileOut.write( '# Time [ms] ' )
            for t in Time:
                FileOut.write( str( t ) + ' ' )
            FileOut.write( '\n' )

        with open( DataFileName, 'a' ) as FileOut:

            # Iterating through an n-dimensional array produces slices along
            # the last axis. This is equivalent to Data[i] in this case

            for TimeSlice in Data:
                FileOut.write( '# New slice\n' )
                np.savetxt( FileOut, TimeSlice )

    return xL, xU, nX, FileArray
