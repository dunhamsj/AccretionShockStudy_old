#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )
from scipy.integrate import trapezoid

from UtilitiesModule import GetFileArray, GetData, ComputeAngleAverage, \
                            Overwrite

#### ========== User Input ==========

rootDirectory = '/lump/data/accretionShockStudy/newData/2D/'

M    = '2.8'
Rpns = '020'
Rs   = '7.00e1'

ID = '2D_M{:}_Rpns{:}_Rs{:}'.format( M, Rpns, Rs )

verbose = False

#### ====== End of User Input =======

Rs = np.float64( Rs )

field = 'LateralMomentumFluxInRadialDirection'

useLogScale = True

ID_GR = 'GR' + ID
ID_NR = 'NR' + ID

# Get mesh
plotFileBaseName = ID_GR + '.plt'
plotFileDirectory = rootDirectory + ID_GR + '/'
plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )
plotFile      = plotFileDirectory + plotFileArray[0]
data, dataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xH, nX \
  = GetData( plotFileDirectory, plotFileBaseName, field+'GR', \
             'spherical', True, argv = [ 'a', plotFileArray[0] ], \
             ReturnTime = False, ReturnMesh = True, Verbose = verbose )

ind = np.where( ( X1[:,0,0] < 0.9 * Rs ) & ( X1[:,0,0] > 0.8 * Rs ) )[0]

OW = Overwrite( 'LatFlux_{:}.dat'.format( ID_GR ) )
if( OW ):

    plotFileBaseNameGR  = ID_GR + '.plt'
    plotFileDirectoryGR = rootDirectory + ID_GR + '/'
    plotFileArrayGR     = GetFileArray \
                            ( plotFileDirectoryGR, plotFileBaseNameGR )
    nSS_GR              = plotFileArrayGR.shape[0]

    plotFileBaseNameNR  = ID_NR + '.plt'
    plotFileDirectoryNR = rootDirectory + ID_NR + '/'
    plotFileArrayNR     = GetFileArray \
                            ( plotFileDirectoryNR, plotFileBaseNameNR )
    nSS_NR              = plotFileArrayNR.shape[0]

    nSS = min( nSS_GR, nSS_NR )

    AA_GR  = np.empty( nSS, np.float64 )
    timeGR = np.empty( nSS, np.float64 )

    AA_NR  = np.empty( nSS, np.float64 )
    timeNR = np.empty( nSS, np.float64 )

    for iSS in range( nSS ):

        if verbose: print('')
        print( '{:}/{:}'.format( iSS+1, nSS ) )

        dataGR, dataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xH, nX, timeGR[iSS] \
          = GetData( plotFileDirectoryGR, plotFileBaseNameGR, field+'GR', \
                     'spherical', True, argv = [ 'a', plotFileArrayGR[iSS] ], \
                     ReturnTime = True, ReturnMesh = True, Verbose = verbose )

        AA = ComputeAngleAverage \
                ( dataGR[ind], X2[0,:,0], dX2[0,:,0] )

        # AA includes \sqrt{\gamma} = \psi^6 * r^2
        AA_GR[iSS] \
          = 4.0 * np.pi \
              * trapezoid \
                  ( AA, \
                    x = X1[ind,0,0], dx = dX1[0,0,0] )

        del dataGR

        dataNR, dataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xH, nX, timeNR[iSS] \
          = GetData( plotFileDirectoryNR, plotFileBaseNameNR, field+'NR', \
                     'spherical', True, argv = [ 'a', plotFileArrayNR[iSS] ], \
                     ReturnTime = True, ReturnMesh = True, Verbose = verbose )

        AA = ComputeAngleAverage \
                ( dataNR[ind], X2[0,:,0], dX2[0,:,0] )

        AA_NR[iSS] \
          = 4.0 * np.pi \
              * trapezoid \
                  ( AA, \
                    x = X1[ind,0,0], dx = dX1[0,0,0] )


        del dataNR

    np.savetxt( 'LatFlux_{:}.dat'.format( ID_GR ), \
                np.vstack( ( timeGR, AA_GR ) ) )
    np.savetxt( 'LatFlux_{:}.dat'.format( ID_NR ), \
                np.vstack( ( timeNR, AA_NR ) ) )

timeGR, dataGR = np.loadtxt( 'LatFlux_{:}.dat'.format( ID_GR ) )
timeNR, dataNR = np.loadtxt( 'LatFlux_{:}.dat'.format( ID_NR ) )

### Plotting

fig, ax  = plt.subplots( 1, 1 )#, figsize = (12,8) )

ax.plot( timeGR[:-1], dataGR[:-1], label = 'GR' )
ax.plot( timeNR[:-1], dataNR[:-1], label = 'NR' )
if( useLogScale ):
    ax.set_yscale( 'symlog', linthresh = 1.0e40 )
    saveFigAs = '/home/kkadoogan/fig.{:}_{:}_SymLog_VsCoordinateTime.png'.format( field, ID )
else:
    saveFigAs = '/home/kkadoogan/fig.{:}_{:}_Linear_VsCoordinateTime.png'.format( field, ID )

ax.legend()

ax.grid()

ax.set_xlabel( r'$\mathrm{Coordinate\ Time}\,\left[\mathrm{ms}\right]$' )

#plt.savefig( saveFigAs, dpi = 300 )
plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
