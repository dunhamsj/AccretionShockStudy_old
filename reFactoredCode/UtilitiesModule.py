#!/usr/bin/env python3

import yt
import numpy as np
from gc import collect

def ComputeAngleAverage( Q, theta, dtheta, dphi ):

    AA = 0.0

    if( Q.ndim > 1 ):

        for iX2 in range( dtheta.shape[0] ):
            for iX3 in range( dphi.shape[0] ):
                AA += Q[iX2,iX3] * np.sin( theta[iX2] ) \
                        * dtheta[iX2] * dphi[iX3]

    else:

        for iX2 in range( dtheta.shape[0] ):
            for iX3 in range( dphi.shape[0] ):
                AA += Q[iX2] * np.sin( theta[iX2] ) \
                        * dtheta[iX2] * dphi[iX3]

    AA  = AA / ( 4.0 * np.pi )

    return AA
# END ComputeAngleAverage


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

        msg = '\n>>>No files found.\n'
        msg += '>>>Double check the path: {:}\n'.format( plotFileDirectory )
        msg += '>>>Double check the PlotFileBaseName: {:}\n' \
               .format( plotFileBaseName )
        msg += '>>>Is it plt_ or just plt?\n'

        assert ( fileArray.shape[0] > 0 ), msg

    return fileArray
# END GetFileArray

def GetData( plotFile, field, verbose = False ):

    """
    GetData
    -------

    str plotFile    : location of plotFile, e.g., /path/to/plotFile.plt########
    str list fields : names of desired fields, e.g., [ 'PF_V1', 'GF_h_2' ]

    """

    yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

    if( verbose ):
        print()
        print( '  Calling GetData...' )
        print( '  ------------------' )
        print( '{:>12} : {:}'.format( 'plotFile', plotFile ) )
        print( '{:>12} : {:}'.format( 'field', field ) )

    ds = yt.load( '{:}'.format( plotFile ) )

    maxLevel_yt = ds.index.max_level
    nX          = ds.domain_dimensions
    xL_yt       = ds.domain_left_edge
    xH_yt       = ds.domain_right_edge

    nDimsX = 1
    if nX[1] > 1: nDimsX += 1
    if nX[2] > 1: nDimsX += 1

    if( nDimsX == 1 ):

        coveringGrid \
          = ds.covering_grid \
              ( level           = maxLevel_yt, \
                left_edge       = xL_yt, \
                dims            = nX * 2**maxLevel_yt, \
                num_ghost_zones = 0 )

    else:

        coveringGrid \
          = ds.covering_grid \
              ( level           = maxLevel_yt, \
                left_edge       = xL_yt, \
                dims            = nX * 2**maxLevel_yt, \
                num_ghost_zones = nX[0] )

        ds.force_periodicity()

    time = np.float64( ds.current_time )

    dataUnits = ''

    if  ( field == 'PF_D' ):

        data = np.copy( coveringGrid['PF_D'].to_ndarray() )
        dataUnits = 'g/cm**3'

    elif( field == 'PF_V1' ):

        data = np.copy( coveringGrid['PF_V1'].to_ndarray() )
        dataUnits = 'km/s'

    elif( field == 'AF_P' ):

        data = np.copy( coveringGrid['AF_P'].to_ndarray() )
        dataUnits = 'erg/cm**3'

    elif( field == 'GF_h_1' ):

        data = np.copy( coveringGrid['GF_h_1'].to_ndarray() )
        dataUnits = ''

    elif( field == 'GF_Gm_11' ):

        data = np.copy( coveringGrid['GF_Gm_11'].to_ndarray() )
        dataUnits = ''

    elif( field == 'GF_Gm_22' ):

        data = np.copy( coveringGrid['GF_Gm_22'].to_ndarray() )
        dataUnits = 'km^2'

    elif( field == 'GF_Alpha' ):

        data = np.copy( coveringGrid['GF_Alph'].to_ndarray() )
        dataUnits = ''

    # --- Derived Fields ---

    elif( field == 'PolytropicConstant' ):

        rho   = np.copy( coveringGrid['PF_D' ].to_ndarray() )
        p     = np.copy( coveringGrid['AF_P' ].to_ndarray() )
        Gamma = np.copy( coveringGrid['AF_Gm'].to_ndarray() )

        data = p / rho**Gamma

        dataUnits = '(erg/cm**3)/(g/cm**3)**(Gamma)'

    elif( field == 'LateralMomentumFluxInRadialDirectionGR' ):

        c = 2.99792458e10

        rho    = np.copy( coveringGrid['PF_D'     ].to_ndarray() )
        e      = np.copy( coveringGrid['PF_E'     ].to_ndarray() )
        p      = np.copy( coveringGrid['AF_P'     ].to_ndarray() )
        V1     = np.copy( coveringGrid['PF_V1'    ].to_ndarray() ) * 1.0e5
        V2     = np.copy( coveringGrid['PF_V2'    ].to_ndarray() )
        Gm11   = np.copy( coveringGrid['GF_Gm_11' ].to_ndarray() )
        Gm22   = np.copy( coveringGrid['GF_Gm_22' ].to_ndarray() ) * (1.0e5)**2
        alpha  = np.copy( coveringGrid['GF_Alpha' ].to_ndarray() )
        SqrtGm = np.copy( coveringGrid['GF_SqrtGm'].to_ndarray() ) * (1.0e5)**2

        h = c**2 + ( e + p ) / rho

        W = 1.0 / np.sqrt( 1.0 - ( Gm11 * V1**2 + Gm22 * V2**2 ) / c**2 )

        data = SqrtGm * rho * Gm22 * V2 * V1 # NR

        data *= alpha * h/c**2 * W**2 # GR corrections

        dataUnits = 'g*cm^2/s^2'

    elif( field == 'LateralMomentumFluxInRadialDirectionNR' ):

        rho    = np.copy( coveringGrid['PF_D'     ].to_ndarray() )
        V1     = np.copy( coveringGrid['PF_V1'    ].to_ndarray() ) * 1.0e5
        V2     = np.copy( coveringGrid['PF_V2'    ].to_ndarray() )
        Gm11   = np.copy( coveringGrid['GF_Gm_11' ].to_ndarray() )
        Gm22   = np.copy( coveringGrid['GF_Gm_22' ].to_ndarray() ) * (1.0e5)**2
        SqrtGm = np.copy( coveringGrid['GF_SqrtGm'].to_ndarray() ) * (1.0e5)**2

        data = SqrtGm * rho * Gm22 * V2 * V1

        dataUnits = 'g*cm^2/s^2'

    else:

        print( '  Invalid field: {:}'.format( field ) )
        print( '  Valid Choices' )
        print( '  -------------' )
        print( '    PF_D' )
        print( '    PF_V1' )
        print( '    AF_P' )
        print( '    GF_h_1' )
        print( '    GF_Gm_11' )
        print( '    GF_Gm_22' )
        print( '    GF_Alpha' )
        print( '    PolytropicConstant' )
        print( '    LateralMomentumFluxInRadialDirectionGR' )
        print( '    LateralMomentumFluxInRadialDirectionNR' )

        exit( '\nExiting...' )

    xL = np.copy( ds.domain_left_edge .to_ndarray() )
    xH = np.copy( ds.domain_right_edge.to_ndarray() )

    dX1 = ( xH[0] - xL[0] ) / np.float64( nX[0] ) * np.ones( nX[0], np.float64 )
    dX2 = ( xH[1] - xL[1] ) / np.float64( nX[1] ) * np.ones( nX[1], np.float64 )
    dX3 = ( xH[2] - xL[2] ) / np.float64( nX[2] ) * np.ones( nX[2], np.float64 )

    X1 = np.linspace( xL[0] + 0.5 * dX1[0], xH[0] - 0.5 * dX1[-1], nX[0] )
    X2 = np.linspace( xL[1] + 0.5 * dX2[0], xH[1] - 0.5 * dX2[-1], nX[1] )
    X3 = np.linspace( xL[2] + 0.5 * dX3[0], xH[2] - 0.5 * dX3[-1], nX[2] )

    if( nDimsX < 3 ):
        dX3[0] = 2.0 * np.pi
        X3 [0] = np.pi

    if( nDimsX < 2 ):
        dX2[0] = np.pi
        X2 [0] = np.pi / 2.0

    # yt has memory leakage issues
    del ds
    collect()

    return time, data, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX
# END GetData

def GetNorm( UseLogScale, Data, vmin = +1.0e100, vmax = -1.0e100, \
             linthresh = 1.0e-2 ):

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, SymLogNorm

    if vmin > +1.0e99: vmin = Data.min()
    if vmax < -1.0e99: vmax = Data.max()

    if UseLogScale:

        if np.any( Data <= 0.0 ):

            Norm = SymLogNorm( vmin = vmin, vmax = vmax, \
                               linthresh = linthresh, base = 10 )

        else:

            Norm = LogNorm   ( vmin = vmin, vmax = vmax )

    else:

        Norm = plt.Normalize ( vmin = vmin, vmax = vmax )

    return Norm
# END GetNorm
