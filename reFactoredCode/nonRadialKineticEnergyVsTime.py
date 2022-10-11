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

field = 'NonRadialKineticEnergyDensity'

useLogScale = True

verbose = False

FileNameGR = 'GR_NonRadKED.dat'
FileNameNR = 'NR_NonRadKED.dat'

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

ind = np.where( ( X1 < 0.90 * Rs ) & ( X1 > 0.80 * Rs ) )[0]

OW = Overwrite( FileNameGR )
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

    E_GR   = np.zeros( nSS_GR, np.float64 )
    timeGR = np.empty( nSS_GR, np.float64 )

    E_NR   = np.zeros( nSS_NR, np.float64 )
    timeNR = np.empty( nSS_NR, np.float64 )

    nSS = min( nSS_GR, nSS_NR )
    for iSS in range( nSS ):

        print( '{:}/{:}'.format( iSS+1, nSS ) )

        plotFileGR = plotFileDirectoryGR + plotFileArrayGR[iSS]
        plotFileNR = plotFileDirectoryNR + plotFileArrayNR[iSS]

        timeGR[iSS], dataGR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFileGR, field+'GR', verbose = verbose )
        t, SqrtGmGR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFileGR, 'GF_SqrtGm', verbose = verbose )

        dX1 *= 1.0e5

        dataGR   = np.copy( dataGR  [ind,:,0] )
        SqrtGmGR = np.copy( SqrtGmGR[ind,:,0] ) * ( 1.0e5 )**2

        t, alpha, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFileGR, 'GF_Alpha', verbose = verbose )
        indd = ( ind[0] + ind[-1] ) // 2
        timeGR[iSS] *= alpha[indd,0,0]

        for iX1 in range( dataGR.shape[0] ):
            for iX2 in range( dataGR.shape[1] ):
                E_GR[iSS] \
                  += 2.0 * np.pi \
                       * dataGR[iX1,iX2] * SqrtGmGR[iX1,iX2] \
                       * dX1[iX1] * dX2[iX2]

        del dataGR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX

        timeNR[iSS], dataNR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFileNR, field+'NR', verbose = verbose )
        t, SqrtGmNR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFileNR, 'GF_SqrtGm', verbose = verbose )

        dX1 *= 1.0e5

        dataNR   = np.copy( dataNR  [ind,:,0] )
        SqrtGmNR = np.copy( SqrtGmNR[ind,:,0] ) * ( 1.0e5 )**2

        for iX1 in range( dataNR.shape[0] ):
            for iX2 in range( dataNR.shape[1] ):
                E_NR[iSS] \
                  += 2.0 * np.pi \
                       * dataNR[iX1,iX2] * SqrtGmNR[iX1,iX2] \
                       * dX1[iX1] * dX2[iX2]

        del dataNR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX

    np.savetxt( FileNameGR, np.vstack( ( timeGR, E_GR ) ) )
    np.savetxt( FileNameNR, np.vstack( ( timeNR, E_NR ) ) )

timeGR, dataGR = np.loadtxt( FileNameGR )
timeNR, dataNR = np.loadtxt( FileNameNR )

### Plotting

fig, ax  = plt.subplots( 1, 1 )#, figsize = (12,8) )

ax.plot( timeGR, dataGR, label = 'GR' )
ax.plot( timeNR, dataNR, label = 'NR' )
if( useLogScale ):
    ax.set_yscale( 'log' )
    saveFigAs = '/home/kkadoogan/fig.{:}_{:}_Log.png'.format( field, ID )
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
