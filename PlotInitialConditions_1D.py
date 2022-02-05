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

Root = '/lump/data/AccretionShockStudy/'

IDs = np.array( [ 'NR1D_M0.14_Mdot0.03_Rs180_PA0.00e-00_nX640', \
                  'GR1D_M0.14_Mdot0.03_Rs180_PA0.00e-00_nX640', \
                  'NR1D_M1.4_Mdot0.3_Rs180_PA0.00e-00_nX640', \
                  'GR1D_M1.4_Mdot0.3_Rs180_PA0.00e-00_nX640' ], str )
c = np.array( [ 'r-', 'b-', 'r--', 'b--' ], str )


Fields = np.array( [ 'PF_D', 'PF_V1', 'AF_P' ], str )

UseLogScales = np.array( [ True, False, True ], bool )

SaveFileAs = 'fig.' + 'IC' + '.png'

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
            axs[iF].plot( r, Data, c[iID], label = ID[0:-23] )
        else:
            axs[iF].plot( r, Data, c[iID] )

        axs[iF].set_ylabel( Field + ' ' + DataUnit )

axs[0].legend()
axs[-1].set_xlabel( r'$r\,\left[\mathrm{km}\right]$' )

plt.savefig( SaveFileAs, dpi = 300 )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
