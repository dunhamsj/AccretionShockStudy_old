#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from UtilitiesModule import GetFileArray, GetData

#### ========== User Input ==========

rootDirectory \
  = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Applications/\
StandingAccretionShock_Relativistic/'

ID = 'GR1D_M3.0_Mdot0.3_Rs060_Rpns020'

field = 'PF_D'

useLogScale = True

saveFigAs = 'fig_{:}.png'.format( ID )

verbose = True

#### ====== End of User Input =======

### Plotting

fig, ax  = plt.subplots( 1, 1, figsize = (12,8) )

plotFileBaseName = ID + '.plt'

plotFileDirectory = rootDirectory + ID + '/'

plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )
plotFile      = plotFileDirectory + plotFileArray[-1]

time, data, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
  = GetData( plotFile, field, verbose = verbose )

ax.plot( X1, data[:,0,0] )

ax.grid()

if( useLogScale ): ax.set_yscale( 'log' )
ax.set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )

#plt.savefig( saveFigAs, dpi = 300 )
plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
