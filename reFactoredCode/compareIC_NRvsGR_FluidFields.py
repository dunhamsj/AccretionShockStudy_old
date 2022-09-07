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
yScales = np.array( [ (2.99792458e10)**(-2), \
                      1.0, \
                      2.99792458e5 ], np.float64 )
useAxis = np.array( [ 0, \
                      0, \
                      1 ], np.int64 )
c = [ [ [ 'r-', 'r--' ], [ 'r.', 'rx' ] ], \
      [ 'b-', 'b--' ] ]

useLogScale = np.array( [ True, False ], bool )

saveFigAs = 'fig_{:}.png'.format( IDs[0] )

verbose = False

#### ====== End of User Input =======

### Plotting

fig = plt.figure( figsize = (12,9) )

me = 10

for i in range( IDs.shape[0] ):

    ax0 = fig.add_subplot( nRows, nCols, i+1 )
    ax1 = ax0.twinx()

    ax0.set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )

    ax0.text( 0.4, 0.7, '{:}'.format( IDs[i][3:] ), fontsize = 15, \
              transform = ax0.transAxes )

    for j in range( fields.shape[0] ):

        plotFileDirectory_NR = rootDirectory_NR
        plotFileDirectory_GR = rootDirectory_GR

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

        if i == 0:

            if( useAxis[j] == 0 ):
                if  ( j == 0 ):
                    ax0.plot( X1, data_NR[:,0,0], c[0][j][0], label = r'$\rho\,c^{2}$ (NR)' )
                    ax0.plot( X1, data_GR[:,0,0], c[0][j][1], label = r'$\rho\,c^{2}$ (GR)' )
                elif( j == 1 ):
                    ax0.plot( X1, data_NR[:,0,0], c[0][j][0], \
                              label = r'$p$ (NR)', markevery = me )
                    ax0.plot( X1, data_GR[:,0,0], c[0][j][1], \
                              label = r'$p$ (GR)', markevery = me )
            elif( useAxis[j] == 1 ):
                ax1.plot( X1, data_NR[:,0,0], c[1][0], label = 'v/c (NR)' )
                ax1.plot( X1, data_GR[:,0,0], c[1][1], label = 'v/c (GR)' )

        else:

            if( useAxis[j] == 0 ):
                if  ( j == 0 ):
                    ax0.plot( X1, data_NR[:,0,0], c[0][j][0] )
                    ax0.plot( X1, data_GR[:,0,0], c[0][j][1] )
                elif( j == 1 ):
                    ax0.plot( X1, data_NR[:,0,0], c[0][j][0], \
                              markevery = me )
                    ax0.plot( X1, data_GR[:,0,0], c[0][j][1], \
                              markevery = me )
            elif( useAxis[j] == 1 ):
                ax1.plot( X1, data_NR[:,0,0], c[1][0] )
                ax1.plot( X1, data_GR[:,0,0], c[1][1] )

    ax0.tick_params( axis = 'x', colors = 'k', labelsize = 15 )
    ax1.tick_params( axis = 'x', colors = 'k', labelsize = 15 )

    ax0.tick_params( axis = 'y', colors = 'r', labelsize = 15 )
    ax1.tick_params( axis = 'y', colors = 'b', labelsize = 15 )

    if( i == 0 ):
        ax0.legend( loc = 3, prop = { 'size' : 15 } )
        ax1.legend( loc = 1, prop = { 'size' : 15 } )

    if( useLogScale[0] ): ax0.set_yscale( 'log' )
    if( useLogScale[1] ): ax1.set_yscale( 'log' )

plt.subplots_adjust( hspace = 0.3 )
#plt.savefig( saveFigAs, dpi = 300 )
plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
