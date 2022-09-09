#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from UtilitiesModule import GetFileArray, GetData

#### ========== User Input ==========

rootDirectory_NR \
  = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Applications/\
StandingAccretionShock_NonRelativistic/'
rootDirectory_GR \
  = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Applications/\
StandingAccretionShock_Relativistic/'

IDs = np.array( [ '1D_M3.0_Mdot0.3_Rs060_Rpns020', \
                  '1D_M3.0_Mdot0.3_Rs075_Rpns020', \
                  '1D_M3.0_Mdot0.3_Rs090_Rpns020' ], str )
R_PNS = 10.0
Rsh = np.array( [ 60.0, 75.0, 90.0 ], np.float64 )

nRows = 2
nCols = 1

yLabels = np.array( [ r'$p/\left(\rho\,c^{2}\right)$', \
                      r'$v/c$' ], str )

saveFigAs = 'fig.CompareNRvsGR_IC_Fluid.png'

verbose = False

#### ====== End of User Input =======

### Plotting

fig = plt.figure( figsize = (12,9) )
fig.suptitle( 'M3.0, Mdot0.3, Rpns010', fontsize = 15, y = 0.92 )

ax0 = fig.add_subplot( nRows, nCols, 1 )
ax1 = fig.add_subplot( nRows, nCols, 2 )

ax0.set_ylabel( yLabels[0] )
ax1.set_ylabel( yLabels[1] )

plotFileDirectory_NR = rootDirectory_NR
plotFileDirectory_GR = rootDirectory_GR

c = [ [ 'r-', 'r--' ], \
      [ 'b-', 'b--' ], \
      [ 'm-', 'm--' ] ]

for i in range( IDs.shape[0] ):

    plotFileBaseName_NR = 'NR' + IDs[i] + '.plt'
    plotFileBaseName_GR = 'GR' + IDs[i] + '.plt'

    plotFileArray_NR \
      = GetFileArray( plotFileDirectory_NR, plotFileBaseName_NR )
    plotFileArray_GR \
      = GetFileArray( plotFileDirectory_GR, plotFileBaseName_GR )

    plotFile_NR = plotFileDirectory_NR + plotFileArray_NR[-1]
    plotFile_GR = plotFileDirectory_GR + plotFileArray_GR[-1]

    time, rho_NR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_NR, 'PF_D', verbose = verbose )
    time, rho_GR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_GR, 'PF_D', verbose = verbose )
    time, p_NR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_NR, 'AF_P', verbose = verbose )
    time, p_GR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_GR, 'AF_P', verbose = verbose )
    time, v_NR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_NR, 'PF_V1', verbose = verbose )
    time, v_GR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_GR, 'PF_V1', verbose = verbose )

    rho_NR *= (2.99792458e10)**2
    rho_GR *= (2.99792458e10)**2

    v_NR /= 2.99792458e5
    v_GR /= 2.99792458e5

    lab = IDs[i][16:21]

    eta = np.copy( ( X1 - R_PNS ) / ( Rsh[i] - R_PNS ) )

    ax0.plot( eta, p_NR[:,0,0] / rho_NR[:,0,0], c[i][0], \
              label = r'\texttt{{NR_'+'{:}}}'.format( lab ) )
    ax0.plot( eta, p_GR[:,0,0] / rho_GR[:,0,0], c[i][1], \
              label = r'\texttt{{GR_'+'{:}}}'.format( lab ) )

    ax1.plot( eta, v_NR[:,0,0], c[i][0] )
    ax1.plot( eta, v_GR[:,0,0], c[i][1] )

ax0.set_xlim( 0.0, 1.0 )
ax1.set_xlim( 0.0, 1.0 )
#ax1.set_xlim( R_PNS, 1.0e2 )

#ax0.set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )
#ax1.set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )
ax0.set_xlabel( r'$\left(r-R_{\textsc{pns}}\right)/\left(R_{s}-R_{\textsc{pns}}\right)$' )
ax1.set_xlabel( r'$\left(r-R_{\textsc{pns}}\right)\left(R_{s}-R_{\textsc{pns}}\right)$' )

ax0.tick_params( axis = 'x', colors = 'k', labelsize = 15 )
ax0.tick_params( axis = 'y', colors = 'k', labelsize = 15 )
ax1.tick_params( axis = 'x', colors = 'k', labelsize = 15 )
ax1.tick_params( axis = 'y', colors = 'k', labelsize = 15 )

ax0.legend( loc = 1, prop = { 'size' : 13 } )

#ax0.set_yscale( 'log' )
ax0.set_ylim( 1.0e-2, 1.5e-1 )
ax1.set_ylim( -6.0e-2, 0.0 )

plt.subplots_adjust( hspace = 0.3 )
plt.savefig( saveFigAs, dpi = 300 )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
