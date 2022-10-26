#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from UtilitiesModule import GetFileArray, GetData

#### ========== User Input ==========

rootDirectory \
  = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Applications/\
StandingAccretionShock_Relativistic/'
#rootDirectory \
#  = '/lump/data/accretionShockStudy/'

ID = '1D_M2.8_Mdot0.3_Rs120'

field = 'AF_Cs'

useLogScale = False

saveFigAs = 'fig.{:}.png'.format( ID )

verbose = True

#### ====== End of User Input =======

### Plotting

fig, ax  = plt.subplots( 1, 1, figsize = (12,8) )

plotFileBaseNameGR = 'GR' + ID + '.plt'
plotFileBaseNameNR = 'NR' + ID + '.plt'

plotFileDirectoryGR = rootDirectory + 'GR' + ID + '/'
plotFileDirectoryNR = rootDirectory + 'NR' + ID + '/'

plotFileArrayGR = GetFileArray( plotFileDirectoryGR, plotFileBaseNameGR )
plotFileArrayNR = GetFileArray( plotFileDirectoryNR, plotFileBaseNameNR )

plotFileGR      = plotFileDirectoryGR + plotFileArrayGR[0]
plotFileNR      = plotFileDirectoryNR + plotFileArrayNR[0]

time, dataGR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
  = GetData( plotFileGR, field, verbose = verbose )
time, dataNR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
  = GetData( plotFileNR, field, verbose = verbose )

ind = np.where( X1 < 115.0 )[0]
print( X1[ind].max() )

ax.plot( X1[ind], dataGR[ind,0,0]/ dataNR[ind,0,0] )

ax.grid()

if( useLogScale ): ax.set_yscale( 'log' )
ax.set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )

#plt.savefig( saveFigAs, dpi = 300 )
plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
