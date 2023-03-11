#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from UtilitiesModule_old import GetData

rootDirectory = '/lump/data/accretionShockStudy/'

# Read in data

R    = np.array( [ 'GR', 'NR' ], str )
M    = np.array( [ '1.4', '2.0', '2.8' ], str )
Rs   = np.array( [ '120', '150', '180' ], str )
Rpns = np.array( [ 4.00e1 ], np.float64 )

R    = np.array( [ 'NR' ], str )
M    = np.array( [ '2.8' ], str )
Rs   = np.array( [ '180' ], str )
Rpns = np.array( [ 4.00e1 ], np.float64 )

arrShape = ( R.shape[0], M.shape[0], Rs.shape[0], Rpns.shape[0] )

v    = np.empty( arrShape, object )
dvdr = np.empty( arrShape, object )
eta  = np.empty( arrShape, object )

for r in range( R.shape[0] ):
    for m in range( M.shape[0] ):
        for rs in range( Rs.shape[0] ):
            for rpns in range( Rpns.shape[0] ):

                ID = '{:}1D_M{:}_Mdot0.3_Rs{:}'.format( R[r], M[m], Rs[rs] )
                plotfileDirectory = rootDirectory + ID + '/'

                plotfileRoot = ID + '.plt_'

                v1, DataUnit, X1_C, X2_C, X3_C, dX1, dX2, dX3, \
                xL, xH, nX, Time \
                  = GetData( plotfileDirectory, plotfileRoot, 'PF_V1', \
                             'spherical', True, MaxLevel = -1, \
                             ReturnTime = True, ReturnMesh = True, \
                             Verbose = True )

                Data, DataUnit, X1_C, X2_C, X3_C, dX1, dX2, dX3, \
                xL, xH, nX, Time \
                  = GetData( plotfileDirectory, plotfileRoot, 'dV1dr', \
                             'spherical', True, MaxLevel = -1, \
                             ReturnTime = True, ReturnMesh = True, \
                             Verbose = True )

                Data = np.copy( Data[1:-1,:,:] )
                v1   = np.copy( v1  [1:-1,:,:] )
                X1_C = np.copy( X1_C[1:-1,:,:] )

                rsh = np.float64( Rs[rs] )
                vsh  = v1[np.where( X1_C < rsh )[0][-1],0,0]

                eta[r,m,rs,rpns] \
                  = ( X1_C[:,0,0] - Rpns[rpns] ) / ( rsh - Rpns[rpns] )

                ind = np.where( eta[r,m,rs,rpns] < 0.99 )[0]

                eta [r,m,rs,rpns] = np.copy( eta[r,m,rs,rpns][ind] )

                Norm = abs( vsh ) / ( rsh - Rpns[rpns] )
                dvdr[r,m,rs,rpns] = np.copy( Data[ind,0,0] ) / Norm
                v   [r,m,rs,rpns] = np.copy( v1  [ind,0,0] )

fig, ax = plt.subplots( 2, 1 )

for r in range( R.shape[0] ):
    for m in range( M.shape[0] ):
        for rs in range( Rs.shape[0] ):
            for rpns in range( Rpns.shape[0] ):
                ax[0].plot( eta[r,m,rs,rpns], v[r,m,rs,rpns] )
                ax[1].plot( eta[r,m,rs,rpns], dvdr[r,m,rs,rpns] )
ax[1].set_xlabel( r'$\eta:=r/R_{\mathrm{S}}$' )
plt.show()
