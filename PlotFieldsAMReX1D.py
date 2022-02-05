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

ID = 'NR1D_M0.14_Mdot0.03_Rs180_PA0.00e-00_nX640'

Field = 'PF_V1'

UseLogScale = False

SaveFileAs = 'fig.' + ID + '.png'

#### ====== End of User Input =======

PlotFileBaseName = ID + '.plt'

DataDirectory += ID + '/'

print( '\nDataDirectory:    {:}'.format( DataDirectory ) )
print( 'PlotFileBaseName: {:}\n'.format( PlotFileBaseName ) )

Data, DataUnit, r, theta, Time, xL, xU \
  = GetData( DataDirectory, PlotFileBaseName, argv, Field, Verbose = True )

### Plotting

Norm = GetNorm( UseLogScale, Data )

fig = plt.figure( figsize = (8,6) )
ax  = fig.add_subplot( 111 )

if( UseLogScale ): ax.set_yscale( 'log' )

ax.plot( r, Data )
ax.set_xlabel( r'$r\,\left[\mathrm{km}\right]$' )

#plt.savefig( SaveFileAs, dpi = 300 )
plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
