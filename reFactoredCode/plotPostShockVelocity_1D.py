#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from UtilitiesModule_old import GetFileArray, GetData

# colorblind-friendly palette: https://gist.github.com/thriveth/8560036
color = ['#377eb8', '#ff7f00', '#4daf4a', \
         '#f781bf', '#a65628', '#984ea3', \
         '#999999', '#e41a1c', '#dede00']

fig, ax  = plt.subplots( 1, 1 )

field = 'PF_V1'
useLogScale = False
saveFigAs = 'fig.PostShockVelocity.pdf'
verbose = True

rootDirectory = '/lump/data/accretionShockStudy/'

M    = [ '1.4'   , '2.8'    ]
Rs   = [ '120'   , '6.00e1' ]
Rpns = [ '4.00e1', '2.00e1' ]

ID_ES = 'GR1D_M{:}_Mdot0.3_Rs{:}'        .format( M[0], Rs[0] )
ID_LS = 'GR1D_M{:}_Mdot0.3_Rs{:}_RPNS{:}'.format( M[1], Rs[1], Rpns[1] )

plotFileBaseNameES = ID_ES + '.plt_'
plotFileBaseNameLS = ID_LS + '.plt'

plotFileDirectoryES = rootDirectory + ID_ES + '/'
plotFileDirectoryLS = rootDirectory + ID_LS + '/'

plotFileArrayES = GetFileArray( plotFileDirectoryES, plotFileBaseNameES )
plotFileArrayLS = GetFileArray( plotFileDirectoryLS, plotFileBaseNameLS )

dataES, DataUnits, X1ES, X2, X3, dX1, dX2, dX3, xL, xH, nX, time \
  = GetData( plotFileDirectoryES, plotFileBaseNameES, field, \
             'spherical', True, argv = [ 'a', '0' ], \
             ReturnTime = True, ReturnMesh = True, Verbose = verbose )

dataLS, DataUnits, X1LS, X2, X3, dX1, dX2, dX3, xL, xH, nX, time \
  = GetData( plotFileDirectoryLS, plotFileBaseNameLS, field, \
             'spherical', True, argv = [ 'a', '0' ], \
             ReturnTime = True, ReturnMesh = True, Verbose = verbose )

etaES \
  = ( X1ES[:,0,0] - np.float64( Rpns[0] ) ) \
  / ( np.float64( Rs[0] ) - np.float64( Rpns[0] ) )

etaLS \
  = ( X1LS[:,0,0] - np.float64( Rpns[1] ) ) \
  / ( np.float64( Rs[1] ) - np.float64( Rpns[1] ) )

indES = np.where( etaES < 0.99 )[0]
indLS = np.where( etaLS < 0.99 )[0]

dataES = np.copy( dataES[indES,0,0] ) / 2.99792458e5
dataLS = np.copy( dataLS[indLS,0,0] ) / 2.99792458e5
etaES  = np.copy( etaES [indES] )
etaLS  = np.copy( etaLS [indLS] )

ax.plot( etaES, dataES, color = color[0], ls = '-', \
        label = r'$\texttt{GR1D\_M1.4\_Rpns040\_Rs120}$' )
ax.plot( etaLS, dataLS, color = color[1], ls = '-', \
        label = r'$\texttt{GR1D\_M2.8\_Rpns020\_Rs060}$' )

ax.legend()
ax.grid()

if( useLogScale ): ax.set_yscale( 'log' )
ax.set_xlabel( r'$\eta:=\left(r-R_{\mathrm{PNS}}\right)/\left(R_{\mathrm{S}}-R_{\mathrm{PNS}}\right)$' )
ax.set_ylabel( r'$v^{1}/c$' )

plt.savefig( saveFigAs, dpi = 300 )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
