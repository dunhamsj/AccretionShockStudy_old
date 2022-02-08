#!/usr/bin/env python3

import numpy as np
import subprocess
from os.path import isfile
from sys import argv, exit
import matplotlib.pyplot as plt
plt.style.use( './Publication.sty' )

from UtilitiesModule import GetData, GetNorm

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True )
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True )
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

#### ========== User Input ==========

#Root = '/lump/data/AccretionShockStudy/'
Root = '/home/dunhamsj/AccretionShockData/'

IDs = np.array( [ 'NR1D_M0.14_Mdot0.03_Rs180_PA0.00e-00_nX640', \
                  'NR1D_M1.4_Mdot0.3_Rs180_PA0.00e-00_nX640' ], str )
c = np.array( [ 'k-', 'k--' ], str )

Fields = np.array( [ 'BernoulliConstant_NR', 'MassConstant_NR' ], str )

UseLogScales = np.array( [ False, True ], bool )

SaveFileAs = 'fig.' + 'Constants' + '.png'

#### ====== End of User Input =======

fig, axs = plt.subplots( Fields.shape[0], 1, figsize = (8,6) )

for iF in range( Fields.shape[0] ):

    Field = Fields[iF]

    for iID in range( IDs.shape[0] ):

        ID = IDs[iID]

        PlotFileBaseName = ID + '.plt'

        DataDirectory = Root + ID + '/'

        Data, DataUnit, r, theta, Time, xL, xU \
          = GetData( DataDirectory, PlotFileBaseName, \
                     argv, Field, Verbose = False )

        Norm = GetNorm( UseLogScales[iF], Data )

        if( UseLogScales[iF] ): axs[iF].set_yscale( 'log' )

        if iF == 0:
            axs[iF].plot( r, Data - Data[0], c[iID], label = ID[0:-23] )
        else:
            axs[iF].plot( r, Data - Data[0], c[iID] )

axs[0].set_ylabel( r'$\mathcal{B}-\mathcal{B}\left(0\right)$' )
axs[0].set_ylabel( r'$(r^{2}\,\rho\,v)-(r^{2}\,\rho\,v)\left(0\right)$' )
axs[0].legend()
axs[-1].set_xlabel( r'$r\,\left[\mathrm{km}\right]$' )

plt.savefig( SaveFileAs, dpi = 300 )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
