#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from UtilitiesModule_old import GetData

rootDirectory = '/lump/data/accretionShockStudy/'

# Read in data

R    = np.array( [ 'NR', 'GR' ], str )
nr   = 0
gr   = 1
M    = np.array( [ '1.4', '2.0', '2.8' ], str )
Rs   = np.array( [ '120', '150', '180' ], str )
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

                Norm = 1.0#abs( vsh ) / ( rsh - Rpns[rpns] )
                dvdr[r,m,rs,rpns] = np.copy( Data[ind,0,0] ) / Norm
                v   [r,m,rs,rpns] = np.copy( v1  [ind,0,0] )

fig, ax = plt.subplots( M.shape[0], Rs.shape[0] )

# colorblind-friendly palette: https://gist.github.com/thriveth/8560036
color = ['#377eb8', '#ff7f00', '#4daf4a', \
         '#f781bf', '#a65628', '#984ea3', \
         '#999999', '#e41a1c', '#dede00']
color = [ 'k' for i in range( 9 ) ]

i = -1
ls = [ '-', '--' ]
for m in range( M.shape[0] ):
    for rs in range( Rs.shape[0] ):
        for rpns in range( Rpns.shape[0] ):
            i += 1
            ax[m,rs].plot( eta[gr,m,rs,rpns], dvdr[gr,m,rs,rpns], \
                           ls = ls[gr], c = color[i], label = R[gr] )
            ax[m,rs].plot( eta[nr,m,rs,rpns], dvdr[nr,m,rs,rpns], \
                           ls = ls[nr], c = color[i], label = R[nr] )
            ax[m,rs].text( 0.4, 0.85, r'$\texttt{{M{:}_Rs{:}}}$' \
                                      .format( M[m], Rs[rs] ), \
                           transform = ax[m,rs].transAxes )
ax[0,0].legend()
fig.supylabel( r'$dv/dr\ \left[\mathrm{s^{-1}}\right]$')# / ( abs(Vsh) / ( Rsh - Rpns ) )$' )
fig.suptitle( 'Suite 1 (Rpns = 40 km)' )
fig.supxlabel( r'$\eta:=\left(r-R_{\mathrm{S}}\right)/\left(R_{\mathrm{PNS}}-R_{\mathrm{S}}\right)$' )
plt.subplots_adjust( wspace = 0.3 )
#plt.show()
plt.savefig( '/home/kkadoogan/fig.velGrad_suite01.png' )
