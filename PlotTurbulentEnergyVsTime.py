#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from MakeDataFile import MakeDataFile

DataDirectory = '/Users/dunhamsj/Research/Data/AccretionShockParameterStudy/'
ID = 'GR2D_M1.4_Mdot0.3_Rs180_PA0.000_nX320x064'
DataDirectory += ID + '/'
PlotFileBaseName = ID + '.plt'

# Get conformal factor

Field = 'GF_Psi'
DataFileName = '.{:}_{:}.dat'.format( ID, Field )
xL, xU, nX, FileArray \
  = MakeDataFile( Field, DataDirectory, DataFileName, PlotFileBaseName )

f = open( DataFileName )
header = f.readline()[16:-2]
FieldUnit = f.readline()[9:-1]
DataShape = ( [ np.int64( dim ) for dim in header.split( ',' ) ] )

Psi \
  = np.loadtxt( DataFileName, dtype = np.float64 ).reshape( DataShape )

# Get turbulent energy density

Field = 'TurbulentEnergyDensity'
DataFileName = '.{:}_{:}.dat'.format( ID, Field )
xL, xU, nX, FileArray \
  = MakeDataFile( Field, DataDirectory, DataFileName, PlotFileBaseName )

f = open( DataFileName )
header = f.readline()[16:-2]
FieldUnit = f.readline()[9:-1]
DataShape = ( [ np.int64( dim ) for dim in header.split( ',' ) ] )

TurbulentEnergyDensity \
  = np.loadtxt( DataFileName, dtype = np.float64 ).reshape( DataShape )

# Get time

Time = list( [ np.float64( t ) for t in f.readline()[12:-2].split(' ') ] )
Time = np.array( Time )

# Get mesh

xL[0] *= 1.0e5
xU[0] *= 1.0e5

dX    = ( xU - xL ) / np.float64( nX )
r     = np.linspace( xL[0] + 0.5 * dX[0], xU[0] - 0.5 * dX[0], nX[0] )
theta = np.linspace( xL[1] + 0.5 * dX[1], xU[1] - 0.5 * dX[1], nX[1] )

# Compute integrand

Integrand = TurbulentEnergyDensity * Psi**6

for i in range( nX[0] ):

    for j in range( nX[1] ):

        Integrand[:,i,j] *= Psi[:,i,j]**6 * r[i]**2 * np.sin( theta[j] )

TurbulentEnergy_X = np.zeros( (DataShape[0],nX[0]), np.float64 )
TurbulentEnergy   = np.zeros( DataShape[0], np.float64 )

# Integrate over angles

print( 'Integrating over angles...' )

for iSS in range( DataShape[0] ):

    for i in range( nX[0] ):

        for j in range( nX[1] - 1 ):

            TurbulentEnergy_X[iSS,i] \
              += 0.5 * ( Integrand[iSS,i,j+1] + Integrand[iSS,i,j] ) * dX[1]

# Integrate over radius

print( 'Integrating over radius...' )

for iSS in range( DataShape[0] ):

    for i in range( nX[0] - 1 ):

        TurbulentEnergy[iSS] \
          += 0.5 * (   TurbulentEnergy_X[iSS,i+1] \
                     + TurbulentEnergy_X[iSS,i  ] ) * dX[0]

TurbulentEnergy *= dX[2]

# Plotting

fig, ax = plt.subplots()
fig.suptitle( ID )
ax.semilogy( Time, TurbulentEnergy )
ax.set_xlabel( 'Time [ms]' )
ax.set_ylabel( 'Turbulent Energy [erg]' )

#plt.show()
plt.savefig( 'fig.TurbulentEnergyVsTime.png', dpi = 300 )

import os
os.system( 'rm -rf __pycache__ ' )
