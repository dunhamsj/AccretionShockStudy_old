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

Rsh = np.array( [ '180', '150', '120', '090', '075', '060' ], str )

IDs = np.array( [ '1D_M2.8_Mdot0.3_Rs{:}_nX380'.format( r ) for r in Rsh ] )
Rsh = np.asarray( Rsh, dtype = np.float64 )

Mpns = '2.8'
Mdot = '0.3'
R_PNS = 10.0

nRows = 3
nCols = 1

yLabels = np.array( [ r'$v/v_{a}$', \
                      r'$\rho/\rho_{a}$', \
                      r'$p/\left(\rho_{a}\,v_{a}^{2}\right)$' ], str )

saveFigAs = 'fig.CompareNRvsGR_IC_KunduCoughlin.png'

verbose = False

#### ====== End of User Input =======

### Plotting

fig, axs = plt.subplots( nRows, nCols, figsize = (12,9) )
fig.suptitle( 'M{:}, Mdot{:}, Rpns{:d}' \
              .format( Mpns, Mdot, np.int64( R_PNS ) ), \
              fontsize = 15, y = 0.92 )

for i in range( nRows ):
    axs[i].set_ylabel( yLabels[i] )
    axs[i].axvline( R_PNS, color = 'k', ls = '--' )

plotFileDirectory_NR = rootDirectory_NR
plotFileDirectory_GR = rootDirectory_GR

c = [ 'red',  \
      'orange',  \
      'yellow',  \
      'green',  \
      'blue',  \
      'magenta' ]

for i in range( IDs.shape[0] ):

    plotFileBaseName_GR = 'GR' + IDs[i] + '.plt'
    plotFileBaseName_NR = 'NR' + IDs[i] + '.plt'

    plotFileArray_GR \
      = GetFileArray( plotFileDirectory_GR, plotFileBaseName_GR )
    plotFileArray_NR \
      = GetFileArray( plotFileDirectory_NR, plotFileBaseName_NR )

    plotFile_GR = plotFileDirectory_GR + plotFileArray_GR[-1]
    plotFile_NR = plotFileDirectory_NR + plotFileArray_NR[-1]

    time, rho_GR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_GR, 'PF_D', verbose = verbose )
    time, rho_NR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_NR, 'PF_D', verbose = verbose )
    time, p_GR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_GR, 'AF_P', verbose = verbose )
    time, p_NR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_NR, 'AF_P', verbose = verbose )
    time, v_GR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_GR, 'PF_V1', verbose = verbose )
    time, v_NR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_NR, 'PF_V1', verbose = verbose )


    lab = IDs[i][16:21]

    #eta = np.copy( ( X1 - R_PNS ) / ( Rsh[i] - R_PNS ) )
    eta = np.copy( X1 )

    #ind = np.where( eta <= R_PNS )[0]
    ind = np.where( eta < Rsh[i] )[0]
    eta = np.copy( eta[ind] )

    indA = np.where( X1 > Rsh[i] )[0][0]

    y0_GR = np.copy( v_GR[ind,0,0] / v_GR[indA,0,0] )
    y0_NR = np.copy( v_NR[ind,0,0] / v_NR[indA,0,0] )

    rhoPrime = rho_GR[indA,0,0] * ( X1 / Rsh[i] )**( -3.0 / 2.0 )

    y1_GR = np.copy( rho_GR[ind,0,0] / rhoPrime )
    y1_NR = np.copy( rho_NR[ind,0,0] / rhoPrime )

    v_GR *= 1.0e5
    v_NR *= 1.0e5

    y2_GR = np.copy( p_GR[ind,0,0] / (rhoPrime*v_GR[indA,0,0]**2) )
    y2_NR = np.copy( p_NR[ind,0,0] / (rhoPrime*v_NR[indA,0,0]**2) )

    sol=2.99792458e10

    axs[0].plot( eta, y0_GR, color = c[i], ls = '-', \
                 label = r'$v_{{a}}={:.2f}$'.format(abs(v_GR[indA,0,0])/sol) )
    axs[0].plot( eta, y0_NR, color = c[i], ls = '--' )
    axs[0].plot( [Rsh[i],Rsh[i]], [0.0,max(y0_GR.max(),y0_NR.max())], \
                 color = c[i] )

    axs[1].plot( eta, y1_GR, color = c[i], ls = '-' )
    axs[1].plot( eta, y1_NR, color = c[i], ls = '--' )
    axs[1].plot( [Rsh[i],Rsh[i]], [0.0,min(y1_GR.min(),y1_NR.min())], \
                 color = c[i] )

    axs[2].plot( eta, y2_GR, color = c[i], ls = '-' )
    axs[2].plot( eta, y2_NR, color = c[i], ls = '--' )
    axs[2].plot( [Rsh[i],Rsh[i]], [0.0,min(y2_GR.min(),y2_NR.min())], \
                 color = c[i] )

for i in range( nRows ):
    axs[i].set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )
#    axs[i].set_xlabel \
#      ( r'$\left(r-R_{\textsc{pns}}\right)/\left(R_{s}-R_{\textsc{pns}}\right)$' )
    axs[i].set_xlim( 7.0, 1000.0 )
    axs[i].set_xscale( 'log' )

    axs[i].tick_params( axis = 'x', colors = 'k', labelsize = 15 )
    axs[i].tick_params( axis = 'y', colors = 'k', labelsize = 15 )

axs[1].set_yscale( 'log' )
axs[2].set_yscale( 'log' )

axs[0].legend( loc = 1, prop = { 'size' : 13 } )

axs[0].set_ylim( -0.01, 0.22 )
axs[1].set_ylim( 1.0e-3, 1.0e6 )
axs[2].set_ylim( 1.0e-5, 1.0e7 )

plt.subplots_adjust( hspace = 0.3 )
#plt.savefig( saveFigAs, dpi = 300 )
plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
