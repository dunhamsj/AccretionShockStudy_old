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

ID = '1D_M1.4_Mdot0.3_Rs180'

field = 'MachNumber'

useLogScale = False

saveFigAs = 'fig.{:}.png'.format( ID )

verbose = True

#### ====== End of User Input =======

### Plotting

fig, ax  = plt.subplots( 1, 1, figsize = (12,8) )

for m in [ '1.4', '2.0', '2.8' ]:
    for rs in [ '120', '150', '180' ]:
        ID = 'GR1D_M{:}_Mdot0.3_Rs{:}'.format( m, rs )
        plotFileBaseNameGR = ID + '.plt'
        plotFileDirectoryGR = rootDirectory + ID + '/'
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
