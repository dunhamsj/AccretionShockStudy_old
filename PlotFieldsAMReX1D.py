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

Root = '/lump/data/AccretionShockStudy/'

IDs = np.array( [ 'GR1D_M1.4_Mdot0.3_Rs150', \
                  'GR1D_M2.0_Mdot0.3_Rs150', \
                  'GR1D_M2.8_Mdot0.3_Rs150', \
                  'NR1D_M1.4_Mdot0.3_Rs150', \
                  'NR1D_M2.0_Mdot0.3_Rs150', \
                  'NR1D_M2.8_Mdot0.3_Rs150' ], str )

Field = 'AF_Cs'

UseLogScale = False

SaveFileAs = 'fig.SoundSpeed.png'

#### ====== End of User Input =======

### Plotting

fig, axs  = plt.subplots( 2, 1, figsize = (12,8) )

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

    ind = np.where( r < 150.0 )[0]

    r    = np.copy( r   [ind] )
    Data = np.copy( Data[ind] )

    Norm = GetNorm( UseLogScale, Data )

    axs[0].plot( r, Data                  , c[i]+ls[i], label = IDs[i] )
    axs[1].plot( r, 2.0 * np.pi * r / Data*1000.0, c[i]+ls[i] )

axs[0].grid()
axs[1].grid()
axs[0].legend()
axs[0].set_ylabel( r'$c_{s}\,\left[\mathrm{km/s}\right]$' )
axs[1].set_ylabel( r'$\tau_{c_{s}}:=\frac{2\pi\,r}{c_{s}}\,\left[\mathrm{ms}\right]$' )
if( UseLogScale ): axs[0].set_yscale( 'log' )
axs[0].set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )
axs[1].set_xlabel( r'Radial Coordinate $\left[\mathrm{km}\right]$' )
plt.subplots_adjust( hspace = 0.3 )

plt.savefig( SaveFileAs, dpi = 300 )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
