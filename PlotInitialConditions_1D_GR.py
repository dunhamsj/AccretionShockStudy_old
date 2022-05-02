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

Root = '/lump/data/AccretionShockStudy/FineGridRuns_1D/'
#Root = HOME + 'AccretionShockData/'

M_lo = 1.0
M_hi = 3.0
nM   = 11
M = np.linspace( M_lo, M_hi, nM )

Rs_lo = 110.0
Rs_hi = 200.0
nRs   = 10
Rs = np.linspace( Rs_lo, Rs_hi, nRs )

IDs = np.empty( nM*nRs, object )

k = -1
for m in range( M.shape[0] ):
    for rs in range( Rs.shape[0] ):
        k += 1
        IDs[k] = 'GR1D_M{:.1f}_Mdot0.3_Rs{:d}'.format \
                   ( M[m], np.int64( Rs[rs] ) )

Fields = np.array( [ 'PF_D', 'PF_V1', 'AF_P', \
                     'SpecificEnthalpy', 'LorentzFactor', 'GF_Alpha' ], str )

UseLogScales = np.array( [ True, False, True, False, False, False ], bool )

SaveFileAs = 'fig.' + 'IC' + '.png'

#### ====== End of User Input =======

fig, axs = plt.subplots( 3, 2, figsize = (14,6) )

k0 = -1
k1 = -1

for iF in range( Fields.shape[0] ):

    Field = Fields[iF]
    print()
    print( Field )

    if iF % 2 == 0: k0 += 1
    else:           k1 += 1

    iID = -1
    for m in range( M.shape[0] ):

        for rs in range( Rs.shape[0] ):

            iID += 1

#            if not ( m  == 0 or m  == 4 or m  == 6 ): continue
#            if not ( rs == 0 or rs == 4 or rs == 6 ): continue
#            if m % 2  == 0: continue
#            if rs % 2 == 0: continue

            print( '  {:}'.format( IDs[iID] ) )
            ID = IDs[iID]

            PlotFileBaseName = ID + '.plt'

            DataDirectory = Root

            Data, DataUnit, r, theta, Time, xL, xU \
              = GetData( DataDirectory, PlotFileBaseName, \
                         Field, argv = [ 'a', '0' ], Verbose = False )

            if iF % 2 == 0:
                axs[k0,0].plot( r, Data )
                if UseLogScales[iF]: axs[k0,0].set_yscale( 'log' )
                axs[k0,0].set_ylabel( Field )
            else:
                axs[k1,1].plot( r, Data )
                if UseLogScales[iF]: axs[k1,1].set_yscale( 'log' )
                axs[k1,1].set_ylabel( Field )
                axs[k1,1].yaxis.set_ticks_position( 'right' )
                axs[k1,1].yaxis.set_label_position( 'right' )

            del Data, DataUnit, r, theta, Time, xL, xU

plt.savefig( SaveFileAs, dpi = 300 )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
