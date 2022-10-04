#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
#plt.style.use( 'publication.sty' )

from UtilitiesModule import GetFileArray, GetData, GetNorm

#### ========== User Input ==========

rootDirectory \
  = '/home/dunhamsj/Work/thornado_GW/SandBox/AMReX/Applications/\
StandingAccretionShock_NonRelativistic/'

ID = 'NR2D_M2.8_Mdot0.3_Rs120_nX768x256'

field = 'PF_E'

useLogScale = True

saveFigAs = 'fig.{:}.png'.format( ID )

verbose = True

#### ====== End of User Input =======

plotFileBaseName = ID + '.plt'

plotFileDirectory = rootDirectory + ID + '/'

plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )
plotFile      = plotFileDirectory + plotFileArray[-1]

time, data, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
  = GetData( plotFile, field, verbose = verbose )
data = np.copy( data[:,:,0] )

Norm = GetNorm( useLogScale, data, vmin = +1.0e100, vmax = -1.0e100, \
                linthresh = 1.0e-2 )

### Plotting

fig = plt.figure( figsize = (16,9) )
ax  = fig.add_subplot( 111, polar = True )
fig.suptitle( ID + '\nTime = ' \
              + '{:.2f} ms'.format( time ) )

im = ax.pcolormesh( X2, X1, data, \
                    cmap = 'viridis', \
                    norm = Norm, \
                    shading = 'nearest' )

ax.set_thetamin( 0.0  )
ax.set_thetamax( 180.0)
ax.set_theta_direction( -1 )
ax.set_theta_zero_location( 'W' )

cbar = fig.colorbar( im )
cbar.set_label( field + ' ' + dataUnits )

plt.savefig( saveFigAs, dpi = 300 )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
