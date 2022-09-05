#!/usr/bin/env python3

import numpy as np

def Overwrite( FileOrDirName, ForceChoice = False, OW = False ):

    if ForceChoice: return OW

    from os.path import isfile, isdir

    OW = True

    if ( isfile( FileOrDirName ) or isdir( FileOrDirName ) ) :

        if ( isdir( FileOrDirName ) and FileOrDirName[-1] != '/' ) :
            FileOrDirName += '/'

        YN = input( '{:} exists. overwrite? (Y/N): '.format( FileOrDirName ) )

        if YN == 'Y' :
            print( 'Overwriting' )
            OW = True
        else:
            print( 'Not overwriting' )
            OW = False

    return OW
# END Overwrite


def GetFileArray( plotFileDirectory, plotFileBaseName ):

    from os import listdir

    fileArray \
      = np.sort( np.array( \
          [ file for file in listdir( plotFileDirectory ) ] ) )

    fileList = []

    for iFile in range( fileArray.shape[0] ):

        sFile = fileArray[iFile]

        if( sFile[0:len(plotFileBaseName)] == plotFileBaseName \
              and sFile[len(plotFileBaseName)+1].isdigit() ) :
            fileList.append( sFile )

    fileArray = np.array( fileList )

    if not fileArray.shape[0] > 0:

        msg = '\n>>>No files found in path {:s}\n'.format( plotFileDirectory )
        msg += '>>>Double check the path.\n'
        msg += '>>>Is it plt_ or just plt?\n'

        assert ( fileArray.shape[0] > 0 ), msg

    return fileArray
# END GetFileArray
