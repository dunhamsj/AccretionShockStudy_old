#!/usr/bin/env python3

import yt
import numpy as np
from os import listdir
from os.path import isfile
from sys import exit

from UtilitiesModule import GetData, GetFileArray, OverwriteFile

def MakeDataFile( Field, DataDirectory, DataFileName, \
                  PlotFileBaseName, UsePhysicalUnits = True, \
                  SSi = -1, SSf = -1, nSS = -1, Verbose = False ):

    if Verbose: print( '\nRunning MakeDataFile...\n' )

    if Verbose: print( 'DataDirectory: {:}\n'.format( DataDirectory ) )

    if( not DataFileName[0] == '.' ): DataFileName = '.' + DataFileName
    if( not DataDirectory[-1] == '/' ): DataDirectory = DataDirectory + '/'

    yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

    c = 2.99792458e10
    if( not UsePhysicalUnits ): c = 1.0

    FileArray = GetFileArray( DataDirectory, PlotFileBaseName, Verbose )

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

    if SSi < 0: SSi = 0
    if SSf < 0: SSf = FileArray.shape[0]
    if nSS < 0: nSS = SSf - SSi

    ind = np.linspace( SSi, SSf, nSS, dtype = np.int64 )

    GenerateFile = True

    if isfile( DataFileName ):

        Overwrite = OverwriteFile( DataFileName )

        GenerateFile = False
        if Overwrite: GenerateFile = True

    if GenerateFile:

        # Put all time-slices into one array to use for movie making

        if   nDimsX == 1: DataShape = (nSS,nX[0])
        elif nDimsX == 2: DataShape = (nSS,nX[0],nX[1])

        Data = np.empty( DataShape, np.float64 )
        Time = np.empty( nSS, np.float64 )

        if Verbose:

            print( 'Generating data file: {:}...'.format( DataFileName ) )

        for i in range( nSS ):

            if i % 10 == 0:
                print( '{:}/{:}'.format( i, nSS ) )

            ds = yt.load( '{:}'.format( DataDirectory + FileArray[ind[i]] ) )

            Data[i], DataUnit, r, theta, t, xL, xU \
              = GetData( DataDirectory, PlotFileBaseName, \
                         [ 'a', FileArray[ind[i]] ], Field )

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
