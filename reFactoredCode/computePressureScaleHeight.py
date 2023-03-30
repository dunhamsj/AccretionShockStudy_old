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
  = '/lump/data/accretionShockStudy/newData/1D/'

ID = 'GR1D_M1.4_Rpns040_Rs1.50e2'

field = 'PressureScaleHeight'

useLogScale = False

saveFigAs = '/home/kkadoogan/fig.{:}.png'.format( ID )

verbose = True

#### ====== End of User Input =======

### Plotting

fig, ax  = plt.subplots( 1, 1 )

plotFileBaseName = ID + '.plt'
plotFileDirectory = rootDirectory + ID + '/'
plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )
plotFile      = plotFileDirectory + plotFileArray[0]
Data, DataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xH, nX, Time \
  = GetData( plotFileDirectory, plotFileBaseName, field, \
             'spherical', True, argv = [ 'a' ], \
             ReturnTime = True, ReturnMesh = True, Verbose = False )

ax.plot( X1[:,0,0], Data[:,0,0] )

ax.grid()

if( useLogScale ): ax.set_yscale( 'log' )
ax.set_title( r'$\texttt{{{:}}}$'.format( ID ) )
ax.set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )
ax.set_ylabel( r'Pressure Scale Height $\left[\mathrm{km}\right]$' )

plt.savefig( saveFigAs, dpi = 300 )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
