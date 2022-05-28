#!/usr/bin/env python3

import numpy as np
import subprocess
from os.path import isfile
from sys import argv, exit
import matplotlib.pyplot as plt
plt.style.use( './Publication.sty' )

from UtilitiesModule import GetData

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True )
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True )
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

#### ========== User Input ==========

#DataDirectory = '/lump/data/AccretionShockStudy/'
#DataDirectory = '/home/dunhamsj/AccretionShockData/'
DataDirectory = '/scratch/dunhamsj/ProductionRuns/'

IDs = np.array( [ 'GR1D_M1.4_Mdot0.3_Rs120', \
                  'GR1D_M1.4_Mdot0.3_Rs150', \
                  'GR1D_M1.4_Mdot0.3_Rs180', \
                  'GR1D_M2.0_Mdot0.3_Rs120', \
                  'GR1D_M2.0_Mdot0.3_Rs150', \
                  'GR1D_M2.0_Mdot0.3_Rs180', \
                  'GR1D_M2.8_Mdot0.3_Rs120', \
                  'GR1D_M2.8_Mdot0.3_Rs150', \
                  'GR1D_M2.8_Mdot0.3_Rs180' ], str )

Field = 'Entropy'

UseLogScale = True

SaveFileAs = 'fig.IC_1D_APS.png'

c = [ 'r-', 'r--', 'r-.', \
      'b-', 'b--', 'b-.', \
      'm-', 'm--', 'm-.' ]

#### ====== End of User Input =======

### Plotting

fig = plt.figure( figsize = (8,6) )
ax  = fig.add_subplot( 111 )

for i in range( IDs.shape[0] ):

    ID = IDs[i]

    PlotFileBaseName = ID + '.plt'

    Root = DataDirectory + ID + '/'

    Data, DataUnit, r, theta, Time, xL, xU \
      = GetData( Root, PlotFileBaseName, argv, Field, Verbose = False )

    lab = ID[5:9] + ID[-6:]

    ax.plot( r, Data, c[i], label = lab )

ax.legend()
if( UseLogScale ): ax.set_yscale( 'log' )

ax.set_xlabel( r'$\mathrm{Radial\ Coordinate}\ \left[\mathrm{km}\right]$' )
ax.set_ylabel( 'Polytropic Constant [cgs]' )

#plt.savefig( SaveFileAs, dpi = 300 )
plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
