#!/usr/bin/env python3

import yt
import numpy as np
from gc import collect

yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

def ComputeTimeScales( plotFileDirectory, rInner, rOuter, rel ):

    # Read in data

    ds       = yt.load( '{:}'.format( plotFileDirectory ) )
    MaxLevel = ds.index.max_level
    nX       = ds.domain_dimensions
    xL       = ds.domain_left_edge
    xU       = ds.domain_right_edge

    CoveringGrid \
      = ds.covering_grid \
          ( level           = MaxLevel, \
            left_edge       = xL, \
            dims            = nX * 2**MaxLevel, \
            num_ghost_zones = nX[0] )

    ds.force_periodicity()

    xL = np.copy( xL.to_ndarray() )
    xU = np.copy( xU.to_ndarray() )

    # Get mesh and isolate elements below shock

    dr = ( xU[0] - xL[0] ) / np.float64( nX[0] )

    r = np.linspace( xL[0] + dr / 2.0, xU[0] - dr / 2.0, nX[0] )

    # Isolate shocked region

    ind = np.where( ( r > rInner ) & ( r < rOuter ) )[0]

    alpha = np.copy( CoveringGrid['GF_Alpha'].to_ndarray()[ind,0,0] )
    Gm11  = np.copy( CoveringGrid['GF_Gm_11'].to_ndarray()[ind,0,0] )
    V     = np.copy( CoveringGrid['PF_V1'   ].to_ndarray()[ind,0,0] )
    Cs    = np.copy( CoveringGrid['AF_Cs'   ].to_ndarray()[ind,0,0] )

    del ds
    collect()

    # Integrate over shocked region to get advection/acoustic times

    tauAd = 0.0
    tauAc = 0.0

    nShocked = np.int64( ind.shape[0] )

    c = 2.99792458e5

    VSq = Gm11 * V  * V

    if  ( rel == 'NR' ):
        lambda0 = V
        lambdap = V + Cs * np.sqrt( 1.0 / Gm11 )
    elif( rel == 'GR' ):
        lambda0 = alpha * V
        lambdap = alpha / ( 1.0 - VSq * Cs**2 / c**4 ) \
                    * ( V * ( 1.0 - Cs**2 / c**2 ) \
                    + Cs * np.sqrt( ( 1.0 - VSq / c**2 ) \
                    * ( 1.0 / Gm11 * ( 1.0 - VSq * Cs**2 / c**4 ) \
                    - V / c * V / c * ( 1.0 - Cs**2 / c**2 ) ) ) )
    for i in range( ind.shape[0] ):

        # Account for integrating inwards in advective timescale
        tauAd += -np.sqrt( Gm11[ind[i]] ) * dr / lambda0[ind[i]]
        tauAc += +np.sqrt( Gm11[ind[i]] ) * dr / lambdap[ind[i]]

    # Convert to ms

    tauAd *= 1.0e3
    tauAc *= 1.0e3
    T_SASI = tauAd + tauAc

    return tauAd, tauAc

if __name__ == '__main__':

    Root = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Applications/'

    rel = 'GR'

    if   rel == 'NR':
        Root += 'StandingAccretionShock_NonRelativistic/'
    elif rel == 'GR':
        Root += 'StandingAccretionShock_Relativistic/'

    M    = np.linspace( 1.4, 2.8, 2, dtype = np.float64 )
    Mdot = np.array( [ 0.3 ], np.float64 )
    Rs   = np.linspace( 30, 180, 64, dtype = np.float64 )
    RPNS = np.linspace( 3 , 42 , 56, dtype = np.float64 )

    for m in M:
        tauAd = []
        tauAc = []
        for mdot in Mdot:
            for rs in Rs:
                tad = []
                tac = []
                for rpns in RPNS:
                    rso = float( rs ) / float( rpns )
                    if rso >= 1.5:
                        ID \
                          = '{:}1D_M{:.1f}_Mdot{:.1f}_Rs{:.3e}_RPNS{:.3e}' \
                            .format( rel, m, mdot, rs, rpns )
                        plotFileDirectory \
                          = Root + '{:}.plt00000000/'.format( ID )

                        rInner = np.float64( rpns )
                        rOuter = np.float64( rs   )

                        tAd, tAc \
                          = ComputeTimeScales \
                              ( plotFileDirectory, rInner, rOuter, rel )
                        print( '{:}, {:.16e}, {:.16e}'.format( ID, tAd, tAc ) )
                        tad.append( tAd )
                        tac.append( tAc )
                    else:
                        tad.append( np.nan )
                        tac.append( np.nan )
                tauAd.append( tad )
                tauAc.append( tac )

        rs = str( [ rs for rs in Rs ] )
        rp = str( [ rpns for rpns in RPNS ] )
        header = '{:}\n{:}'.format( rs, rp )
        tauAd = np.array( tauAd, np.float64 )
        tauAc = np.array( tauAc, np.float64 )
        np.savetxt( 'tauAd_{:}_M{:.1f}.dat' \
                    .format( rel, m ), tauAd, header = header )
        np.savetxt( 'tauAc_{:}_M{:.1f}.dat' \
                    .format( rel, m ), tauAc, header = header )

    import os
    os.system( 'rm -rf __pycache__ ' )
