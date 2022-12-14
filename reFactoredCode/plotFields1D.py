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
  = '/lump/data/accretionShockStudy/'

ID = '1D_M2.8_Mdot0.3_Rs9.00e1_RPNS2.00e1'

field = 'PF_E'

useLogScale = True

saveFigAs = 'fig.{:}.png'.format( ID )

verbose = True

#### ====== End of User Input =======

### Plotting

fig, ax  = plt.subplots( 1, 1, figsize = (12,8) )

plotFileBaseNameGR = 'GR' + ID + '.plt'

plotFileDirectoryGR = rootDirectory + 'GR' + ID + '_nX0512/'

plotFileArrayGR = GetFileArray( plotFileDirectoryGR, plotFileBaseNameGR )
plotFileGR      = plotFileDirectoryGR + plotFileArrayGR[0]

time, dataGR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
  = GetData( plotFileGR, field, verbose = verbose )

ax.plot( X1, dataGR[:,0,0] )

ax.grid()

if( useLogScale ): ax.set_yscale( 'log' )
ax.set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )

#plt.savefig( saveFigAs, dpi = 300 )
plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
