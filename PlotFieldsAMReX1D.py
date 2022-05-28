#!/usr/bin/env python3

import numpy as np
import subprocess
from os.path import isfile
from sys import argv, exit
import matplotlib.pyplot as plt

from UtilitiesModule import GetData, GetNorm

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True )
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True )
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

#### ========== User Input ==========

Root = '/lump/data/AccretionShockStudy/'

IDs = np.array( [ 'GR1D_M2.0_Mdot0.1_Rs150', \
                  'GR1D_M2.0_Mdot0.3_Rs150', \
                  'GR1D_M2.0_Mdot0.5_Rs150' ], str )

Field = 'AF_Cs'

UseLogScale = False

SaveFileAs = 'fig.' + IDs[0] + '.png'

#### ====== End of User Input =======

### Plotting

fig = plt.figure( figsize = (8,6) )
ax  = fig.add_subplot( 111 )

for i in range( IDs.shape[0] ):

    PlotFileBaseName = IDs[i] + '.plt'

    DataDirectory = Root + IDs[i] + '/'

    print( '\nDataDirectory:    {:}'.format( DataDirectory ) )
    print( 'PlotFileBaseName: {:}\n'.format( PlotFileBaseName ) )

    Data , DataUnit, r, theta, phi, dr, dtheta, dphi, xL, xU, nX, Time \
      = GetData( DataDirectory, PlotFileBaseName, Field, 'spherical', True,
                 argv = [ 'a', '00000000' ], Verbose = True, \
                 ReturnTime = True, ReturnMesh = True )

    Norm = GetNorm( UseLogScale, Data )

    ax.plot( r, Data, label = IDs[i] )

ax.grid()
ax.legend()
if( UseLogScale ): ax.set_yscale( 'log' )
ax.set_xlabel( r'$r\,\left[\mathrm{km}\right]$' )

#plt.savefig( SaveFileAs, dpi = 300 )
plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
