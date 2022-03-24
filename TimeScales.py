#!/usr/bin/env python3

import yt
import numpy as np

yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

class TimeScales:

    def __init__( self ):
        return

    def ComputeTimeScales( self, DataDirectory, rInner, rOuter ):

        # Read in data

        ds       = yt.load( '{:}'.format( DataDirectory ) )
        MaxLevel = ds.index.max_level
        nX       = ds.domain_dimensions
        xL       = ds.domain_left_edge
        xU       = ds.domain_right_edge

        CoveringGrid \
          = ds.covering_grid \
              ( level           = MaxLevel, \
                left_edge       = xL, \
                dims            = nX * 2**MaxLevel, \
                num_ghost_zones = 0 )

        xL = xL.to_ndarray()
        xU = xU.to_ndarray()

        # Get mesh and isolate elements below shock

        dr = ( xU[0] - xL[0] ) / np.float64( nX[0] )

        r = np.linspace( xL[0] + dr / 2.0, xU[0] - dr / 2.0, nX[0] )

        # Isolate shocked region

        ind = np.where( ( r > rInner ) & ( r < rOuter ) )[0]

        V1 = CoveringGrid['PF_V1'].to_ndarray()[ind,0,0]
        Cs = CoveringGrid['AF_Cs'].to_ndarray()[ind,0,0]

        # Integrate over shocked region to get advection/acoustic times

        tauAd = 0.0
        tauAc = 0.0

        nShocked = np.int64( ind.shape[0] )

        for i in range( ind.shape[0] ):

            tauAd += dr / np.abs( V1[ind[i]] )
            tauAc += dr / ( Cs[ind[i]] - np.abs( V1[ind[i]] ) )

        # Convert to ms

        tauAd *= 1.0e3
        tauAc *= 1.0e3
        T_SASI = tauAd + tauAc

        return T_SASI

if __name__ == '__main__':

    Root = '/lump/data/AccretionShockStudy/'

    M    = np.array( [ '1.4', '2.0', '2.8' ], str )
    Mdot = np.array( [ '0.3' ], str )
    Rs   = np.array( [ '120', '150', '180' ], str )

    TS = TimeScales()

    for m in range( M.shape[0] ):
        for mdot in range( Mdot.shape[0] ):
            for rs in range( Rs.shape[0] ):

                ID \
                  = 'GR1D_M{:}_Mdot{:}_Rs{:}'.format( M[m], Mdot[mdot], Rs[rs] )
                DataDirectory = Root + '{:}/'.format( ID )
                DataDirectory += '{:}.plt_00000000/'.format( ID )

                rInner = 4.00e1
                rOuter = np.float64( Rs[rs] )

                print( '\nM{:}_Mdot{:}_Rs{:}'.format \
                       ( M[m], Mdot[mdot], Rs[rs] ) )
                T_SASI = TS.ComputeTimeScales( DataDirectory, rInner, rOuter )
                print( 'T_SASI: {:.3e}'.format( T_SASI ) )

    import os
    os.system( 'rm -rf __pycache__ ' )
