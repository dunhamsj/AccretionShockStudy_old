#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from UtilitiesModule import GetFileArray, GetData

#### ========== User Input ==========

#rootDirectory \
#  = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Applications/\
#StandingAccretionShock_Relativistic/'
rootDirectory \
  = '/lump/data/accretionShockStudy/newData/resolutionStudy_earlyStage/'

ID = 'GR1D_M1.4_Rpns040_Rs1.20e2'
#ID = 'GR1D_M2.8_Rpns020_Rs6.00e1'

field = 'PressureScaleHeight'

useLogScale = False

saveFigAs = '/home/kkadoogan/fig.{:}.png'.format( ID )

verbose = True

#### ====== End of User Input =======

### Plotting

fig, ax  = plt.subplots( 1, 1 )

dr = [ 1.0, 0.5, 0.25, 0.125 ]
#dr = [ 0.5, 0.25, 0.125, 0.0625 ]
i = -1
for nx in [ '0140', '0280', '0560', '1120' ]:

    IDD = ID + '_nX{:}'.format( nx )

    plotFileDirectory = rootDirectory + IDD + '/'
    plotFileBaseName = IDD + '.plt'
    plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )
    plotFile      = plotFileDirectory + plotFileArray[0]
    Data, DataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xH, nX, Time \
      = GetData( plotFileDirectory, plotFileBaseName, field, \
                 'spherical', True, argv = [ 'a' ], \
                 ReturnTime = True, ReturnMesh = True, Verbose = False )

    ind = np.where( X1[:,0,0] < 1.19e2 )[0]
    ind = np.copy( ind[1:] )
    i += 1
    ax.plot( X1[ind,0,0], Data[ind,0,0], \
             label = r'$dr={:}\,\mathrm{{km}}$' \
                     .format( dr[i] ) )

ax.grid()

ax.legend()
if( useLogScale ): ax.set_yscale( 'log' )
ax.set_title( r'$\texttt{{{:}}}$'.format( ID ) )
ax.set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )
#ax.set_ylabel( r'Pressure Scale Height $\left[\mathrm{km}\right]$' )
ax.set_ylabel( r'$p\left|\frac{dp}{dr}\right|^{-1}\,\left[\mathrm{km}\right]$' )

plt.savefig( saveFigAs, dpi = 300 )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
