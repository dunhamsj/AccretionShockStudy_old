#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

from UtilitiesModule import GetNorm, GetFileArray
from MakeDataFile import MakeDataFile, ReadHeader

ID = 'GR1D_M2.0_Mdot0.3_Rs150.nX0640'

plotFileDirectory = '/lump/data/accretionShockStudy/'
plotFileDirectory += ID + '/'

plotFileBaseName = ID + '.plt'

field = 'PolytropicConstant'

useLogScale = True

SSi = -1 # -1 -> SSi = 0
SSf = -1 # -1 -> plotFileArray.shape[0] - 1
nSS = 100 # -1 -> plotFileArray.shape[0]

maxLevel = -1

showIC = True

plotMesh = False

verbose = True

movieRunTime = 10.0 # seconds

useCustomLimits = False

#### ====== End of User Input =======

dataFileDirectory = '.{:s}_{:s}_MovieData1D'.format( ID, field )
movieName         = 'mov.{:s}_{:s}.mp4'.format( ID, field )

# Append "/" if not present
if( not plotFileDirectory[-1] == '/' ): plotFileDirectory += '/'
if( not dataFileDirectory[-1] == '/' ): dataFileDirectory += '/'

MakeDataFile( field, plotFileDirectory, dataFileDirectory, \
              plotFileBaseName, \
              SSi = SSi, SSf = SSf, nSS = nSS, \
              maxLevel = maxLevel, verbose = verbose, \
              forceChoice = False, OW = True )

# Ignore last file (t=0 without perturbation)
plotFileArray \
  = GetFileArray \
      ( plotFileDirectory, plotFileBaseName, \
        SSi = SSi, SSf = SSf, nSS = nSS )#[:-1]

if SSi < 0: SSi = 0
if SSf < 0: SSf = plotFileArray.shape[0] - 1
if nSS < 0: nSS = plotFileArray.shape[0]

fig = plt.figure()
ax  = fig.add_subplot( 111 )

if useCustomLimits: ax.set_ylim( ymin, ymax )
ymin = +np.inf
ymax = -np.inf
for j in range( nSS ):
    iSS = SSi + np.int64( ( SSf - SSi ) / ( nSS - 1 ) * j )
    dataFile = dataFileDirectory + plotFileArray[iSS] + '.dat'
    data = np.loadtxt( dataFile )
    ymin = min( ymin, data.min() )
    ymax = max( ymax, data.max() )
if not useCustomLimits: ax.set_ylim( ymin, ymax )

def f(t):

    iSS = SSi + np.int64( ( SSf - SSi + 1 ) / nSS ) * t

    dataFile = dataFileDirectory + plotFileArray[iSS] + '.dat'

    dataShape, dataUnits, Time, X1_C, X2_C, X3_C, dX1, dX2, dX3 \
      = ReadHeader( dataFile )

    data = np.loadtxt( dataFile ).reshape( np.int64( dataShape ) )

    return data, dataUnits, X1_C, dX1, Time

data0, dataUnits, X1_C0, dX10, Time0 = f(0)

xL = np.array( [ X1_C0[0 ]-0.5*dX10[0 ] ], np.float64 )
xU = np.array( [ X1_C0[-1]+0.5*dX10[-1] ], np.float64 )

Norm = GetNorm( useLogScale, data0, vmin = ymin, vmax = ymax )

TimeUnits   = 'ms'
LengthUnits = 'km'

ax.set_xlim( xL[0], xU[0] )

ax.set_xlabel( 'Isotropic Radius [km]' )
ax.set_ylabel( field + ' ' + dataUnits )

time_text = plt.text( 0.5, 0.7, '', transform = ax.transAxes )

if( useLogScale ): ax.set_yscale( 'log' )

ms = 5
line, = ax.plot( [],[], 'k.', markersize = ms, label = 'u(t)' )
if showIC: IC, = ax.plot( [],[], 'r.', markersize = ms, label = 'u(0)' )
if plotMesh: mesh, = ax.plot( [],[], 'b.', label = 'mesh' )

def InitializeFrame():

    line.set_data([],[])
    time_text.set_text('')
    if showIC: IC.set_data([],[])
    if plotMesh: mesh.set_data([],[])

    if showIC and plotMesh: ret = ( line, time_text, IC, mesh )
    elif showIC:            ret = ( line, time_text, IC )
    elif plotMesh:          ret = ( line, time_text, mesh )
    else:                   ret = ( line, time_text )

    return ret

def UpdateFrame( t ):

    data, dataUnits, X1_C, dX1, Time = f(t)

    line.set_data( X1_C, data.flatten() )
    time_text.set_text( 'Time = {:.3e} {:}'.format( Time, TimeUnits ) )
    if showIC: IC.set_data( X1_C0, data0.flatten() )
    if plotMesh: mesh.set_data( X1_C - 0.5 * dX1, 0.5 * ( ymin + ymax ) )

    if showIC and plotMesh: ret = ( line, time_text, IC, mesh )
    elif showIC:            ret = ( line, time_text, IC )
    elif plotMesh:          ret = ( line, time_text, mesh )
    else:                   ret = ( line, time_text )

    return ret

ax.legend()

anim = animation.FuncAnimation( fig, UpdateFrame, \
                                init_func = InitializeFrame, \
                                frames = nSS, \
                                blit = True )

fps = max( 1, nSS / movieRunTime )

print( 'Saving movie {:}...'.format( movieName ) )

anim.save( movieName, fps = fps, dpi = 300 )

print( 'Done!' )

import os
os.system( 'rm -rf __pycache__ ' )
