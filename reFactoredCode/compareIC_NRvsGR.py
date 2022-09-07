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
nRows = 3
nCols = 1

fields = np.array( [ 'PF_D', \
                     'AF_P', \
                     'PF_V1' ], str )
yLabels = np.array( [ r'$\rho\,\left[\mathrm{g\,cm}^{-3}\right]$', \
                      r'$P\,\left[\mathrm{erg\,cm}^{-3}\right]$', \
                      r'$v/c$' ], str )
yScales = np.array( [ 1.0, \
                      1.0, \
                      2.99792458e5 ], np.float64 )
useLogScale = np.array( [ True, \
                          True, \
                          False ], bool )

saveFigAs = 'fig.CompareNRvsGR_IC_Fluid.png'

verbose = False

#### ====== End of User Input =======

### Plotting

fig = plt.figure( figsize = (12,9) )
fig.suptitle( 'M3.0, Mdot0.3, Rpns020', fontsize = 15 )

c = [ [ 'r-', 'r--' ], \
      [ 'b-', 'b--' ], \
      [ 'm-', 'm--' ] ]

for j in range( fields.shape[0] ):

    ax0 = fig.add_subplot( nRows, nCols, j+1 )

    ax0.set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )
    ax0.set_ylabel( yLabels[j] )

    ax0.set_xlim( 0.0, 1.0e2 )

    plotFileDirectory_NR = rootDirectory_NR
    plotFileDirectory_GR = rootDirectory_GR

    for i in range( IDs.shape[0] ):

#        ax0.text( 0.4, 0.7, '{:}'.format( IDs[i][3:] ), fontsize = 15, \
#                  transform = ax0.transAxes )

        plotFileBaseName_NR = 'NR' + IDs[i] + '.plt'
        plotFileBaseName_GR = 'GR' + IDs[i] + '.plt'

        plotFileArray_NR \
          = GetFileArray( plotFileDirectory_NR, plotFileBaseName_NR )
        plotFileArray_GR \
          = GetFileArray( plotFileDirectory_GR, plotFileBaseName_GR )

        plotFile_NR = plotFileDirectory_NR + plotFileArray_NR[-1]
        plotFile_GR = plotFileDirectory_GR + plotFileArray_GR[-1]

        time, data_NR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFile_NR, fields[j], verbose = verbose )
        time, data_GR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFile_GR, fields[j], verbose = verbose )

        data_NR /= yScales[j]
        data_GR /= yScales[j]

        if( j == 2 ):
            lab = IDs[i][16:21]
            ax0.plot( X1, data_NR[:,0,0], c[i][0], label = 'NR_'+'{:}'.format( lab ) )
            ax0.plot( X1, data_GR[:,0,0], c[i][1], label = 'GR_'+'{:}'.format( lab ) )
        else:
            ax0.plot( X1, data_NR[:,0,0], c[i][0] )
            ax0.plot( X1, data_GR[:,0,0], c[i][1] )

    ax0.tick_params( axis = 'x', colors = 'k', labelsize = 15 )
    ax0.tick_params( axis = 'y', colors = 'k', labelsize = 15 )

    if( j == 2 ):
        ax0.legend( loc = 2, prop = { 'size' : 13 } )

    if( useLogScale[j] ): ax0.set_yscale( 'log' )

plt.subplots_adjust( hspace = 0.3 )
#plt.savefig( saveFigAs, dpi = 300 )
plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
