#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use( 'publication.sty' )

from UtilitiesModule import GetNorm, MapCenterToCorners
from MakeDataFile import MakeDataFile, ReadHeader

"""

Creates a directory with structure as laid out
in MakeDataFile.py and makes a movie from it

Usage:
  $ python3 makeMovie2D.py

"""

#### ========== User Input ==========

nPanels = 2

ID                = np.empty( nPanels, object )
plotfileDirectory = np.empty( nPanels, object )
plotfileBaseName  = np.empty( nPanels, object )
field             = np.empty( nPanels, object )
cmap              = np.empty( nPanels, object )
useCustomLimits   = np.empty( nPanels, bool )
vmin              = np.empty( nPanels, np.float64 )
vmax              = np.empty( nPanels, np.float64 )
dataDirectory     = np.empty( nPanels, object )
plotfileArray     = np.empty( nPanels, object )
dataUnits         = np.empty( nPanels, object )
time              = np.empty( nPanels, object )
X1_C              = np.empty( nPanels, object )
X2_C              = np.empty( nPanels, object )
dX1               = np.empty( nPanels, object )
dX2               = np.empty( nPanels, object )
data              = np.empty( nPanels, object )
x1L               = np.empty( nPanels, np.float64 )
x1H               = np.empty( nPanels, np.float64 )
x2L               = np.empty( nPanels, np.float64 )
x2H               = np.empty( nPanels, np.float64 )
X1c               = np.empty( nPanels, object )
X2c               = np.empty( nPanels, object )
norm              = np.empty( nPanels, object )
im                = np.empty( nPanels, object )
cbar              = np.empty( nPanels, object )

# ID to be used for naming purposes
ID[0] = 'NR2D_M1.4_Rpns040_Rs150_Mdot0.3'
ID[1] = 'NR2D_M1.4_Rpns040_Rs180_Mdot0.3'

title = 'NR2D_M1.4_Rpns040_Mdot0.3'

# Directory containing AMReX plotfiles
plotfileDirectory[0] = '/lump/data/accretionShockStudy/\
newRuns/newProductionRuns/{:}/'.format( ID[0] )
plotfileDirectory[1] = '/lump/data/accretionShockStudy/\
newRuns/newProductionRuns/{:}/'.format( ID[1] )

# plotfile base name (e.g., Advection2D.plt######## -> Advection2D.plt )
for i in range( nPanels ):
    plotfileBaseName[i] = ID[i] + '.plt'

# Field to plot
field[0] = 'PolytropicConstant'
field[1] = 'PolytropicConstant'

# Plot data in log10-scale?
useLogScale = True

# Only use every <plotEvery> plotfile
plotEvery = 10

# Colormap
cmap[0] = 'RdBu'
cmap[1] = 'RdBu'

# First and last snapshots and number of snapshots to include in movie
SSi = -1 # -1 -> SSi = 0
SSf = -1 # -1 -> plotfileArray.shape[0] - 1
nSS = -1 # -1 -> plotfileArray.shape[0]

Verbose = True

for i in range( nPanels ):
    useCustomLimits[i] = True
    vmin[i] = -1.0e-6
    vmax[i] = +1.0e-6

MovieRunTime = 10.0 # seconds

#### ====== End of User Input =======

polar = True

for i in range( nPanels ):
    dataDirectory[i] = '.{:}/'.format( ID[i] )

MovieName = 'mov.{:}_{:}_{:}_{:}.mp4' \
            .format( ID[0], field[0], ID[1], field[1] )

TimeUnits = 'ms'
X1Units   = 'km'
X2Units   = 'rad'

nSSS = 100000
for i in range( nPanels ):
    plotfileArray[i] \
      = MakeDataFile( field[i], plotfileDirectory[i], dataDirectory[i], \
                      plotfileBaseName[i], 'spherical', \
                      SSi = SSi, SSf = SSf, nSS = nSS, \
                      UsePhysicalUnits = True, \
                      MaxLevel = -1, Verbose = Verbose )
    plotfileArray[i] = np.copy( plotfileArray[i][:-1] )
    plotfileArray[i] = np.copy( plotfileArray[i][::plotEvery] )

    nSSS = min( nSSS, plotfileArray[i].shape[0] )

if nSS < 0:
    nSS = nSSS

# Ensure all arrays have same number of plotfiles
for i in range( nPanels ):
    plotfileArray[i] = np.copy( plotfileArray[i][0:nSS] )

def f(t):

    for i in range( nPanels ):

        fileDirectory = dataDirectory[i] + plotfileArray[i][t] + '/'

        TimeFile = fileDirectory + '{:}.dat'.format( 'Time' )
        X1File   = fileDirectory + '{:}.dat'.format( 'X1' )
        X2File   = fileDirectory + '{:}.dat'.format( 'X2' )
        dX1File  = fileDirectory + '{:}.dat'.format( 'dX1' )
        dX2File  = fileDirectory + '{:}.dat'.format( 'dX2' )
        DataFile = fileDirectory + '{:}.dat'.format( field[i] )

        DataShape, dataUnits[i], MinVal, MaxVal = ReadHeader( DataFile )

        time[i] = np.loadtxt( TimeFile )
        X1_C[i] = np.loadtxt( X1File   ).reshape( np.int64( DataShape ) )
        X2_C[i] = np.loadtxt( X2File   ).reshape( np.int64( DataShape ) )
        dX1 [i] = np.loadtxt( dX1File  ).reshape( np.int64( DataShape ) )
        dX2 [i] = np.loadtxt( dX2File  ).reshape( np.int64( DataShape ) )
        data[i] = np.loadtxt( DataFile )

    return np.copy( data ), dataUnits, X1_C, X2_C, dX1, dX2, np.copy( time )

Data0, DataUnits, X1_C, X2_C, dX1, dX2, Time = f(0)
Data1, DataUnits, X1_C, X2_C, dX1, dX2, Time = f(0)

for i in range( nPanels ):

    x1L[i] = X1_C[i][0 ,0 ] - 0.5 * dX1[i][0 ,0 ]
    x1H[i] = X1_C[i][-1,-1] + 0.5 * dX1[i][-1,-1]
    x2L[i] = X2_C[i][0 ,0 ] - 0.5 * dX2[i][0 ,0 ]
    x2H[i] = X2_C[i][-1,-1] + 0.5 * dX2[i][-1,-1]

    if not useCustomLimits[i]:
        vmin[i] = +np.inf
        vmax[i] = -np.inf
        for j in range( nSS ):
            DataFile \
              = dataDirectory[i] + plotfileArray[i][j] \
                  + '/{:}.dat'.format( field[i] )
            DataShape, DataUnits, MinVal, MaxVal = ReadHeader( DataFile )
            vmin[i] = min( vmin[i], MinVal )
            vmax[i] = max( vmax[i], MaxVal )

fig = plt.figure()
ax  = fig.add_subplot( 111, polar = polar )
fig.suptitle( r'$\texttt{{{:}}}$'.format( title ) )

for i in range( nPanels ):

    X1c[i], X2c[i] = MapCenterToCorners( X1_C[i], X2_C[i], dX1[i], dX2[i] )

    if i == 1:

        X1c[i] = np.copy( X1c[i][:,0] )
        X2c[i] = 2.0 * np.pi - np.copy( X2c[i][0,:] )

    else:

        X1c[i] = np.copy( X1c[i][:,0] )
        X2c[i] = np.copy( X2c[i][0,:] )

    X1c[i], X2c[i] = np.meshgrid( X2c[i], X1c[i] )

ax.grid( False )

ax.set_thetamin( 0.0 )
ax.set_thetamax( 360.0 )
ax.set_rmin( 0.0 )
rmax  = 0.0
for i in range( nPanels ):
    rmax = max( x1H[i], rmax )
ax.set_rmax( rmax )

ax.set_theta_zero_location( 'N' ) # z-axis vertical
ax.set_theta_direction( -1 )

if field[0] == field[1] and not useCustomLimits[0]:

    vmn = np.min( vmin )
    vmx = np.max( vmax )

    mn = vmn
    mx = vmx

    if( ( vmx >= 0.0 ) & ( vmn < 0.0 ) ):

        mn = min( vmn, -vmx )
        mx = max( vmx, -vmn )

    for i in range( nPanels ):
        vmin[i] = mn
        vmax[i] = mx

time_text = ax.set_title( '' )

for i in range( nPanels ):
    d = ( Data1[i] - Data0[i] ) / Data0[i] + 1.0e-17
    d[0] -= 2.0e-17

    norm[i] = GetNorm( useLogScale, d, vmin = vmin[i], vmax = vmax[i] )

    im[i] = ax.pcolormesh( X1c[i], X2c[i], d, \
                           cmap = cmap[i], \
                           norm = norm[i], \
                           shading = 'flat' )

    if i == 0:
        location = 'right'
        labelpad = 0
    else:
        location = 'left'
        labelpad = 0

    cbar[i] = fig.colorbar( im[i], location = location )
    cbar[i].set_label( r'$dK/K$' + ' {:}'.format( ID[i][-13:-8] ) )
      #cbar[i].set_label \
      #  ( field[i] + ' ' + r'$\mathrm{{{:}}}$' \
      #                    .format( dataUnits[i][1:-1] ), labelpad = labelpad )

def InitializeFrame():

    ret = []
    Data, DataUnits, X1_C, X2_C, dX1, dX2, Time = f(0)
    for i in range( nPanels ):
        Data[i] = ( Data[i] - Data0[i] ) / Data0[i] + 1.0e-17
        im[i].set_array( Data[i].flatten() )
        ret.append( im[i] )

    time_text.set_text('')
    ret.append( time_text )

    return ret

def UpdateFrame(t):

    print( '    {:}/{:}'.format( t, nSS ) )

    ret = []
    Data, DataUnits, X1_C, X2_C, dX1, dX2, Time = f(t)
    for i in range( nPanels ):
        Data[i] = ( Data[i] - Data0[i] ) / Data0[i] + 1.0e-17
        im[i].set_array( Data[i].flatten() )
        ret.append( im[i] )

    time_text.set_text( r'$t={:.3e}\ \left[\mathrm{{{:}}}\right]$' \
                        .format( Time[i], TimeUnits ) )
    ret.append( time_text )

    return ret

# Call the animator
anim \
  = animation.FuncAnimation \
      ( fig, UpdateFrame, init_func = InitializeFrame, \
        frames = nSS, blit = True)

fps = max( 1, nSS // MovieRunTime )

print( '\n  Making Movie' )
print( '  ------------' )

anim.save( MovieName, fps = fps, dpi = 300 )

import os
os.system( 'rm -rf __pycache__ ' )
