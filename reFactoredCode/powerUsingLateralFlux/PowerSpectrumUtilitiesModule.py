#!/usr/bin/env python3

from gc import collect
from psutil import Process
from os import getpid
import yt
import numpy as np
from scipy.integrate import simps

def printMemUsage( verbose = True ):

    process = Process( getpid() )
    memory = process.memory_info().rss / 1024.0

    if verbose:
        print( 'mem: {:.3e} kB'.format( memory ) )
        return
    else:
        return memory
# END printMemUsage

def ReadFields( plotFile, field, verbose = False ):

    """
    ReadFields
    ----------

    str plotFile    : location of plotFile, e.g., /path/to/plotFile.plt########
    str list fields : names of desired fields, e.g., [ 'PF_V1', 'GF_h_2' ]

    """

    yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

    if( verbose ):
        print()
        print( '  Calling ReadFields...' )
        print( '  ---------------------' )
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

    data = np.empty( (2,nX[0],nX[1],nX[2]), np.float64 )

    data[0] = np.copy( coveringGrid['GF_Psi'].to_ndarray() )

    if( field == 'DivV2' ):

        data[1] = np.copy( coveringGrid['PF_V2'].to_ndarray() )

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

        d = SqrtGm * rho * Gm22 * V2 * V1 # NR

        data[1] = d * alpha * h/c**2 * W**2 # GR corrections

    elif( field == 'LateralMomentumFluxInRadialDirectionNR' ):

        rho    = np.copy( coveringGrid['PF_D'     ].to_ndarray() )
        V1     = np.copy( coveringGrid['PF_V1'    ].to_ndarray() ) * 1.0e5
        V2     = np.copy( coveringGrid['PF_V2'    ].to_ndarray() )
        Gm11   = np.copy( coveringGrid['GF_Gm_11' ].to_ndarray() )
        Gm22   = np.copy( coveringGrid['GF_Gm_22' ].to_ndarray() ) * (1.0e5)**2
        SqrtGm = np.copy( coveringGrid['GF_SqrtGm'].to_ndarray() ) * (1.0e5)**2

        data[1] = SqrtGm * rho * Gm22 * V2 * V1

    else:

        print( 'PowerSpectrumUtilitiesModule' )
        print( '----------------------------' )
        print( 'Invalid choice of field: {:}\n'.format( field ) )
        print( 'Valid choices' )
        print( '-------------' )
        print( '  DivV2' )
        print( '  LateralMomentumFluxInRadialDirectionGR' )
        print( '  LateralMomentumFluxInRadialDirectionNR' )
        exit( 'Exiting...' )

    xL = ds.domain_left_edge .to_ndarray()
    xH = ds.domain_right_edge.to_ndarray()

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

    return time, data, X1, X2, X3, dX1, dX2, dX3, nX
# END ReadFields


def GetMesh( plotFile, verbose = False ):

    """
    GetMesh
    -------

    str plotFile : location of plotFile, e.g., /path/to/plotFile.plt########

    """

    yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

    if( verbose ):
        print()
        print( '  Calling ReadField...' )
        print( '  --------------------' )
        print( '{:>12} : {:}'.format( 'plotFile', plotFile ) )

    ds = yt.load( '{:}'.format( plotFile ) )

    nX = ds.domain_dimensions
    xL = ds.domain_left_edge .to_ndarray()
    xH = ds.domain_right_edge.to_ndarray()

    dX1 = ( xH[0] - xL[0] ) / np.float64( nX[0] )
    dX2 = ( xH[1] - xL[1] ) / np.float64( nX[1] )
    dX3 = ( xH[2] - xL[2] ) / np.float64( nX[2] )

    X1 = np.linspace( xL[0] + 0.5 * dX1, xH[0] - 0.5 * dX1, nX[0] )
    X2 = np.linspace( xL[1] + 0.5 * dX2, xH[1] - 0.5 * dX2, nX[1] )
    X3 = np.linspace( xL[2] + 0.5 * dX3, xH[2] - 0.5 * dX3, nX[2] )

    # yt has memory leakage issues
    del ds
    collect()

    return X1, X2, X3, dX1, dX2, dX3, nX
# END GetMesh

if __name__ == '__main__':

    # Test ReadFields

    plotFile \
      = '/lump/data/accretionShockStudy/\
NR2D_M2.0_Mdot0.3_Rs120/NR2D_M2.0_Mdot0.3_Rs120.plt_00412668'
    field = [ 'PF_V1', 'PF_D' ]
    start = printMemUsage( verbose = False )
    for i in range( 1 ):
        t, F, X1, X2, X3, dX1, dX2, dX3, nX \
          = ReadFields( plotFile, field, verbose = True )
    stop = printMemUsage( verbose = False )
    print( stop - start )

    # Test GetMesh

    plotFile \
      = '/lump/data/accretionShockStudy/\
NR2D_M2.0_Mdot0.3_Rs120/NR2D_M2.0_Mdot0.3_Rs120.plt_00412668'
    start = printMemUsage( verbose = False )
    for i in range( 1 ):
        X1, X2, X3, dX1, dX2, dX3, nX = GetMesh( plotFile, verbose = True )
    print( X1 )
    print( X2 )
    print( X3 )
    print( dX1 )
    print( dX2 )
    print( dX3 )
    stop = printMemUsage( verbose = False )
    print( stop - start )
