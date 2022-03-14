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

IDs = np.array( [ 'GR1D_M0.14_Mdot0.03_Rs180_PA0.00e-00_nX640', \
                  'GR1D_M0.14_Mdot0.03_Rs180_PA0.00e-00_nX640' ], str )

c = np.array( [ 'rs', 'b.', 'r.-', 'b.-' ], str )

Fields = np.array( [ 'Entropy', 'PF_V1', 'AF_P' ], str )

UseLogScales = np.array( [ True, False, True ], bool )

SaveFileAs = 'fig.' + 'IC' + '.png'

#### ====== End of User Input =======

fig, axs = plt.subplots( Fields.shape[0], 2, figsize = (14,6) )

Data = np.empty( (IDs.shape[0],Fields.shape[0],640), np.float64 )
DataUnit = []

for iF in range( Fields.shape[0] ):

    Field = Fields[iF]

    for iID in range( IDs.shape[0] ):

        ID = IDs[iID]

        PlotFileBaseName = ID + '.plt'

        suffix = ''
        if iID == 1: suffix = '_old'
        DataDirectory = Root + ID + suffix + '/'

        Data[iID,iF], DataUnitt, r, theta, Time, xL, xU \
          = GetData( DataDirectory, PlotFileBaseName, \
                     argv, Field, Verbose = False )

        if iID == 0: DataUnit.append( DataUnitt )

DataUnit = np.array( DataUnit )

for iF in range( Fields.shape[0] ):

    Field = Fields[iF]

    for iID in range( IDs.shape[0] ):

        ID = IDs[iID]

        Norm = GetNorm( UseLogScales[iF], Data[iID] )

        if( UseLogScales[iF] ): axs[iF,0].set_yscale( 'log' )

        if iF == 0:
            if iID == 0:
                axs[iF,0].plot( r, Data[iID,iF], c[iID], markerfacecolor = 'none', label = 'new' )
            else:
                axs[iF,0].plot( r, Data[iID,iF], c[iID], label = 'old' )
        elif iF == 1:
            axs[iF,0].plot( r, Data[iID,iF], c[iID] )
        else:
            axs[iF,0].plot( r, Data[iID,iF], c[iID] )

    axs[0,0].set_ylabel( 'Polytropic Constant [cgs]' )
    axs[1,0].set_ylabel( 'V [km/s]' )
    axs[2,0].set_ylabel( 'P [erg/cc]' )

    axs[iF,1].semilogy( r, np.abs( ( Data[0,iF] - Data[1,iF] ) / Data[0,iF] ) )

axs[0,1].set_ylabel( r'$|(K_{N}-K_{O})/K_{N}|$' )
axs[1,1].set_ylabel( r'$|(v_{N}-v_{O})/v_{N}|$' )
axs[2,1].set_ylabel( r'$|(p_{N}-p_{O})/p_{N}|$' )

axs[0,0].legend()
axs[-1,0].set_xlabel( r'$r\,\left[\mathrm{km}\right]$' )
axs[-1,1].set_xlabel( r'$r\,\left[\mathrm{km}\right]$' )

plt.savefig( SaveFileAs, dpi = 300 )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
