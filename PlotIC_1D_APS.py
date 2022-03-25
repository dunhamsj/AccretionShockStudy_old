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

DataDirectory = '/lump/data/AccretionShockStudy/'
#DataDirectory = '/home/dunhamsj/AccretionShockData/'
#DataDirectory = '/scratch/dunhamsj/ProductionRuns/'

ID = 'GR1D_M1.4_Mdot0.3_Rs120'

Fields = np.array( [ 'PF_D', 'PF_V1', 'AF_P' ], str )

UseLogScale = np.array( [ True, False, True ], bool )

SaveFileAs = 'fig.IC_1D_APS.png'

#### ====== End of User Input =======

### Plotting

PlotFileBaseName = ID + '.plt'

Root = DataDirectory + ID + '/'

D, DUnit, r, theta, Time, xL, xU \
  = GetData( Root, PlotFileBaseName, argv, Fields[0], Verbose = False )
V, VUnit, r, theta, Time, xL, xU \
  = GetData( Root, PlotFileBaseName, argv, Fields[1], Verbose = False )
P, PUnit, r, theta, Time, xL, xU \
  = GetData( Root, PlotFileBaseName, argv, Fields[2], Verbose = False )

fig, ax1 = plt.subplots( 1, 1, figsize = (8,6) )
ax2 = ax1.twinx()

c = 2.99792458e10
l1D, = ax1.plot( r, D*c**2, 'r-', label = r'$\rho\,c^{2}$' )
l1P, = ax1.plot( r, P, 'r--', label = r'$p$' )
l2,  = ax2.plot( r, V, 'b-' )

ax2.spines['left'].set_color('red')
ax2.spines['right'].set_color('blue')

size = 20

ax1.tick_params( axis = 'y', colors = 'red'  , labelsize = size )
ax2.tick_params( axis = 'y', colors = 'blue' , labelsize = size )
ax1.tick_params( axis = 'x', colors = 'black', labelsize = size )

ax1.set_yscale( 'log' )

plt.legend( [l1D,l1P,l2], [r'$\rho\,c^{2}$',r'$p$',r'$v$'], \
            prop = {'size':15} )

ax1.set_xlabel( r'$\mathrm{Radial\ Coordinate}\ \left[\mathrm{km}\right]$', fontsize = size )
ax1.set_ylabel( r'$\mathrm{erg\,cm}^{-3}$', color = 'red' , fontsize = size )
ax2.set_ylabel( r'$\mathrm{km\,s}^{-1}$'  , color = 'blue', fontsize = size )

plt.savefig( SaveFileAs, dpi = 300, bbox_inches = 'tight' )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
