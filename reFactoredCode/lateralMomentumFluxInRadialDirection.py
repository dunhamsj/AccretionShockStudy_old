#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from UtilitiesModule import GetFileArray, GetData, ComputeAngleAverage, \
                            Overwrite

#### ========== User Input ==========

rootDirectory \
  = '/lump/data/accretionShockStudy/'
#rootDirectory \
#  = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Applications/\
#StandingAccretionShock_Relativistic/'

Rs = 120

ID = '2D_M2.8_Mdot0.3_Rs{:d}'.format( Rs )

field = 'LateralMomentumFluxInRadialDirection'

useLogScale = True

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
  = GetData( plotFile, field+'GR', verbose = True )

ind = np.where( X1 < 0.95 * Rs )[0][-1]

OW = Overwrite( 'GR_LatFlux.dat' )
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
          = GetData( plotFileGR, field+'GR', verbose = verbose )
        t, alphaGR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFileGR, 'GF_Alpha', verbose = verbose )
        timeGR[iSS] *= alphaGR[ind,0,0]
        AA_GR[iSS] = ComputeAngleAverage( dataGR[ind], X2, dX2, dX3 )

        del dataGR

        timeNR[iSS], dataNR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFileNR, field+'NR', verbose = verbose )
        AA_NR[iSS] = ComputeAngleAverage( dataNR[ind], X2, dX2, dX3 )

        del dataNR

    np.savetxt( 'GR_LatFlux.dat', np.vstack( ( timeGR, AA_GR ) ) )
    np.savetxt( 'NR_LatFlux.dat', np.vstack( ( timeNR, AA_NR ) ) )

timeGR, dataGR = np.loadtxt( 'GR_LatFlux.dat' )
timeNR, dataNR = np.loadtxt( 'NR_LatFlux.dat' )

### Plotting

fig, ax  = plt.subplots( 1, 1 )#, figsize = (12,8) )

ax.plot( timeGR, dataGR, label = 'GR' )
ax.plot( timeNR, dataNR, label = 'NR' )
if( useLogScale ):
    ax.set_yscale( 'symlog', linthresh = 1.0e40 )
    saveFigAs = '/home/kkadoogan/fig.{:}_{:}_SymLog.png'.format( field, ID )
else:
    saveFigAs = '/home/kkadoogan/fig.{:}_{:}_Linear.png'.format( field, ID )

ax.legend()

ax.grid()

ax.set_xlabel( r'$\mathrm{Proper\ Time}\,\left[\mathrm{ms}\right]$' )

plt.savefig( saveFigAs, dpi = 300 )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
