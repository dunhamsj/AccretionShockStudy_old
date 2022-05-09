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

DataDirectory = '/lump/data/AccretionShockStudy/'

ID = 'GR1D_M2.0_Mdot0.3_Rs150_entropyPert_PA1.00e-06'

Field = 'PF_D'

UseLogScale = False

SaveFileAs = 'fig.' + ID + '.png'

#### ====== End of User Input =======

PlotFileBaseName = ID + '.plt'

DataDirectory += ID + '/'

print( '\nDataDirectory:    {:}'.format( DataDirectory ) )
print( 'PlotFileBaseName: {:}\n'.format( PlotFileBaseName ) )

Data , DataUnit, r, theta, phi, dr, dtheta, dphi, xL, xU, nX, Time \
  = GetData( DataDirectory, PlotFileBaseName, Field, 'spherical', True,
             argv = [ 'a', '00000000' ], Verbose = True, \
             ReturnTime = True, ReturnMesh = True )
Data0, DataUnit, r, theta, phi, dr, dtheta, dphi, xL, xU, nX, Time \
  = GetData( DataDirectory, PlotFileBaseName, Field, 'spherical', True, \
             argv = [ 'a', '99999999'], Verbose = True, \
             ReturnTime = True, ReturnMesh = True )

Data = ( Data - Data0 ) / Data0

### Plotting

Norm = GetNorm( UseLogScale, Data )

fig = plt.figure( figsize = (8,6) )
ax  = fig.add_subplot( 111 )

if( UseLogScale ): ax.set_yscale( 'log' )

ax.plot( r, Data, 'k.' )
ax.set_xlabel( r'$r\,\left[\mathrm{km}\right]$' )

#plt.savefig( SaveFileAs, dpi = 300 )
plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
