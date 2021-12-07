#!/usr/bin/env python3

import yt
import numpy as np

yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

def ComputeTimeScales( DataDirectory, rInner, rOuter ):

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
        tauAc += dr / np.abs( Cs[ind[i]] )

    # Convert to ms

    tauAd *= 1.0e3
    tauAc *= 1.0e3

    print( 'Range of integration: [ {:.2e}, {:.2e} ] km'.format \
             ( rInner, rOuter ) )

    print( 'tauAd = {:.3e} ms'.format( tauAd ) )
    print( 'tauAc = {:.3e} ms'.format( tauAc ) )

    return

ID = 'GR1D_M1.4_Mdot0.3_Rs180_PA0.000_nX640'

DataDirectory = "/Users/dunhamsj/Research/Data/"
DataDirectory += "AccretionShockParameterStudy/{:}/".format( ID )
DataDirectory += "{:}.plt_00000000/".format( ID )

rInner = 4.00e1
rOuter = 1.80e2

ComputeTimeScales( DataDirectory, rInner, rOuter )

import os
os.system( 'rm -rf __pycache__ ' )
