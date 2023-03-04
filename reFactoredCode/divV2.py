#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )
from scipy.integrate import trapezoid

from UtilitiesModule import GetFileArray, GetData, ComputeAngleAverage, \
                            Overwrite

#### ========== User Input ==========

rootDirectory \
  = '/lump/data/accretionShockStudy/'
#rootDirectory \
#  = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Applications/\
#StandingAccretionShock_Relativistic/'

Rs = '120'

ID = '2D_M2.8_Mdot0.3_Rs{:}'.format( Rs )

Rs = np.float64( Rs )

field = 'DivV2'

useLogScale = False

verbose = False

#### ====== End of User Input =======

ID_GR = 'GR' + ID
ID_NR = 'NR' + ID

# Get mesh
plotFileBaseName = ID_GR + '.plt'
plotFileDirectory = rootDirectory + ID_GR + '/'
plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )
plotFile      = plotFileDirectory + plotFileArray[0]
time, data, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
  = GetData( plotFile, field, verbose = True )

ind = np.where( ( X1 < 0.9 * Rs ) & ( X1 > 0.8 * Rs ) )[0]

OW = Overwrite( 'DivV2_{:}.dat'.format( ID_GR ) )
if( OW ):

    plotFileBaseNameGR  = ID_GR + '.plt'
    plotFileDirectoryGR = rootDirectory + ID_GR + '/'
    plotFileArrayGR     = GetFileArray( plotFileDirectoryGR, plotFileBaseNameGR )
    nSS_GR              = plotFileArrayGR.shape[0]

    plotFileBaseNameNR  = ID_NR + '.plt'
    plotFileDirectoryNR = rootDirectory + ID_NR + '/'
    plotFileArrayNR     = GetFileArray( plotFileDirectoryNR, plotFileBaseNameNR )
    nSS_NR              = plotFileArrayNR.shape[0]

    AA_GR  = np.empty( nSS_GR, np.float64 )
    timeGR = np.empty( nSS_GR, np.float64 )

    AA_NR  = np.empty( nSS_NR, np.float64 )
    timeNR = np.empty( nSS_NR, np.float64 )

    nSS = min( nSS_GR, nSS_NR )
    for iSS in range( nSS ):

        print( '{:}/{:}'.format( iSS+1, nSS ) )

        plotFileGR = plotFileDirectoryGR + plotFileArrayGR[iSS]
        plotFileNR = plotFileDirectoryNR + plotFileArrayNR[iSS]

        timeGR[iSS], dataGR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFileGR, field, verbose = verbose )
        AA = ComputeAngleAverage \
               ( dataGR[ind], X2, dX2, dX3, \
                 nX = [ dataGR[ind].shape[0], \
                        dataGR[ind].shape[1], \
                        dataGR[ind].shape[2] ] )
        AA_GR[iSS] \
          = 4.0 * np.pi \
              * trapezoid \
                  ( AA * X1[ind]**2, \
                    x = X1[ind], dx = dX1[0] )

        del dataGR

        timeNR[iSS], dataNR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFileNR, field, verbose = verbose )
        AA = ComputeAngleAverage \
               ( dataNR[ind], X2, dX2, dX3, \
                 nX = [ dataNR[ind].shape[0], \
                        dataNR[ind].shape[1], \
                        dataNR[ind].shape[2] ] )
        AA_NR[iSS] \
          = 4.0 * np.pi \
              * trapezoid \
                  ( AA * X1[ind]**2, \
                    x = X1[ind], dx = dX1[0] )

        del dataNR

    np.savetxt( 'DivV2_{:}.dat'.format( ID_GR ), \
                np.vstack( ( timeGR, AA_GR ) ) )
    np.savetxt( 'DivV2_{:}.dat'.format( ID_NR ), \
                np.vstack( ( timeNR, AA_NR ) ) )

timeGR, dataGR = np.loadtxt( 'DivV2_{:}.dat'.format( ID_GR ) )
timeNR, dataNR = np.loadtxt( 'DivV2_{:}.dat'.format( ID_NR ) )

### Plotting

fig, ax  = plt.subplots( 1, 1 )#, figsize = (12,8) )

ax.plot( timeGR, dataGR, label = 'GR' )
ax.plot( timeNR, dataNR, label = 'NR' )
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
