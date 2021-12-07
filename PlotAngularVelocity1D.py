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

DataDirectory = HOME + 'Research/Data/AccretionShockParameterStudy/'

ID = 'GR1D_M1.4_Mdot0.3_Rs180_PA0.000_nX320'

Field = 'PF_V2'

SaveFileAs = 'fig.' + ID + '.png'

#### ====== End of User Input =======

PlotFileBaseName = ID + '.plt'

DataDirectory += ID + '_Nodal/'

print( '\nDataDirectory:    {:}'.format( DataDirectory ) )
print( 'PlotFileBaseName: {:}\n'.format( PlotFileBaseName ) )

Data, DataUnit, r, theta, Time, xL, xU \
  = GetData( DataDirectory, PlotFileBaseName, argv, Field )

### Plotting

fig = plt.figure( figsize = (8,6) )
ax  = fig.add_subplot( 111 )

ind = np.where( ( r[:,0] < 1.15e2 ) & ( r[:,0] > 1.05e2 ) )[0]

for i in range( ind.shape[0] ):
    if i % 1 == 0:
        ax.plot( theta[ind[i]], Data[ind[i]], label = 'r = {:.3e} km'.format( r[ind[i],0] ) )
ax.set_ylabel( r'$v^2\,\left[\mathrm{rad/s}\right]$' )
ax.legend()
ax.set_xlabel( r'$\theta$' )

#plt.savefig( SaveFileAs, dpi = 300 )
plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
