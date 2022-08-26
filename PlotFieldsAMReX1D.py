#!/usr/bin/env python3

import numpy as np
import subprocess
from os.path import isfile
from sys import argv, exit
import matplotlib.pyplot as plt
plt.style.use( 'Publication.sty' )

from UtilitiesModule import GetData, GetNorm

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True )
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True )
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

#### ========== User Input ==========

#Root = '/lump/data/AccretionShockStudy/'
#Root = THORNADO_DIR + 'SandBox/AMReX/Euler_Relativistic_IDEAL/'
Root = '/home/kkadoogan/'

IDs = np.array( [ 'NR1D_M2.0_Mdot0.3_Rs120' ], str )

Field = 'PF_D'

UseLogScale = True

SaveFileAs = '/home/kkadoogan/fig.png'

#### ====== End of User Input =======

### Plotting

fig, ax  = plt.subplots( 1, 1, figsize = (12,8) )

ls = [ '-', '-', '-', '--', '--', '--' ]
c = [ 'r', 'b', 'm', 'r', 'b', 'm' ]

for i in range( IDs.shape[0] ):

    PlotFileBaseName = IDs[i] + '.plt'

    DataDirectory = Root + IDs[i] + '/'

    print( '\nDataDirectory:    {:}'.format( DataDirectory ) )
    print( 'PlotFileBaseName: {:}\n'.format( PlotFileBaseName ) )

    Data , DataUnit, r, theta, phi, dr, dtheta, dphi, xL, xU, nX, Time \
      = GetData( DataDirectory, PlotFileBaseName, Field, 'spherical', True,
                 argv = [ 'a', '00000000' ], Verbose = True, \
                 ReturnTime = True, ReturnMesh = True )

#    ind = np.where( r < 150.0 )[0]
#
#    r    = np.copy( r   [ind] )
#    Data = np.copy( Data[ind] )

    Norm = GetNorm( UseLogScale, Data )

    ax.plot( r, Data, c[i]+ls[i], label = IDs[i] )

ax.grid()
if( UseLogScale ): ax.set_yscale( 'log' )
ax.set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )

#plt.savefig( SaveFileAs, dpi = 300 )
plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
