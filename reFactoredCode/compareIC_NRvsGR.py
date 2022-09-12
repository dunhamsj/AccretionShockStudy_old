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

#IDs = np.array( [ '1D_M2.8_Mdot0.3_Rs060_nX380', \
#                  '1D_M2.8_Mdot0.3_Rs075_nX380', \
#                  '1D_M2.8_Mdot0.3_Rs090_nX380', \
#                  '1D_M2.8_Mdot0.3_Rs120_nX380', \
#                  '1D_M2.8_Mdot0.3_Rs150_nX380', \
#                  '1D_M2.8_Mdot0.3_Rs180_nX380' ], str )

IDs = np.array( [ '1D_M3.0_Mdot0.3_Rs020_nX400' ], str )

R_PNS = 10.0
#Rsh = np.array( [ 60.0, 75.0, 90.0, 120.0, 150.0, 180.0 ], np.float64 )
Rsh = np.array( [ 20.0 ], np.float64 )

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

c = [ 'red',  \
      'orange',  \
      'yellow',  \
      'green',  \
      'blue',  \
      'magenta' ]

y0Min = +np.inf
y0Max = -np.inf
y1Min = +np.inf
y1Max = -np.inf
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
    #eta = np.copy( X1 )

    ind = np.where( eta <= 1.0 )[0]
    eta = np.copy( eta[ind] )

    y0_NR = np.copy( p_NR[ind,0,0] / rho_NR[ind,0,0] )
    y0_GR = np.copy( p_GR[ind,0,0] / rho_GR[ind,0,0] )
    y1_NR = np.copy( v_NR[ind,0,0] )
    y1_GR = np.copy( v_GR[ind,0,0] )

    ax0.plot( eta, y0_NR, color = c[i], ls = '-', \
              label = r'\texttt{{NR_'+'{:}}}'.format( lab ) )
    ax0.plot( eta, y0_GR, color = c[i], ls = '--', \
              label = r'\texttt{{GR_'+'{:}}}'.format( lab ) )

    ax1.plot( eta, y1_NR, color = c[i], ls = '-' )
    ax1.plot( eta, y1_GR, color = c[i], ls = '--' )

#ax0.set_xlim( R_PNS, 5.0e1 )
#ax1.set_xlim( R_PNS, 5.0e1 )

#ax0.set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )
#ax1.set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )
ax0.set_xlabel \
  ( r'$\left(r-R_{\textsc{pns}}\right)/\left(R_{s}-R_{\textsc{pns}}\right)$' )
ax1.set_xlabel \
  ( r'$\left(r-R_{\textsc{pns}}\right)\left(R_{s}-R_{\textsc{pns}}\right)$' )

ax0.tick_params( axis = 'x', colors = 'k', labelsize = 15 )
ax0.tick_params( axis = 'y', colors = 'k', labelsize = 15 )
ax1.tick_params( axis = 'x', colors = 'k', labelsize = 15 )
ax1.tick_params( axis = 'y', colors = 'k', labelsize = 15 )

ax0.legend( loc = 1, prop = { 'size' : 13 } )

plt.subplots_adjust( hspace = 0.3 )
plt.savefig( saveFigAs, dpi = 300 )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
