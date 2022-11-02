#!/usr/bin/env python3

import numpy as np
import subprocess
from os.path import isfile
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use( 'publication.sty' )

from UtilitiesModule import GetFileArray
from MakeDataFile import MakeDataFile, ReadHeader

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True )
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True )
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

#### ========== User Input ==========

# ID to be used for naming purposes
ID = 'GR1D_M2.8_Mdot0.3_Rs9.00e1_RPNS2.00e1'

useTwoAxes = True

tb = 0.0#2.503e2

# Directory containing AMReX plotfiles
plotFileDirectory0 = '/lump/data/accretionShockStudy/'
plotFileDirectory0 += ID + '_nX0512/'

plotFileDirectory1 = plotFileDirectory0

# PlotFile base name (e.g., Advection2D.plt######## -> Advection2D.plt )
plotFileBaseName = ID + '.plt'

def setFields( plotWhat ):

    if plotWhat == 'hydro':

        # Hydro

        useTwoAxes = True

        # Fields to plot
        fields0 = np.array( [ 'PF_D', 'AF_P' ] )
        yScale0 = np.array( [ (2.99792458e10)**(-2), 1.0 ], np.float64 )
        label0 = [ r'$\rho\,c^{2}$', r'$p$' ]
        useCustomY0Limits = False
        y0Min = -1.0e23
        y0Max = 1.0e1
        useLogScale0 = True
        y0Label = r'$\mathrm{erg\,cm}^{-3}$'

        fields1 = np.array( [ 'PF_V1' ] )
        yScale1 = np.array( [ 2.99792458e5 ], np.float64 )
        label1 = [ r'$v^{1}/c$' ]
        useCustomY1Limits = False
        y1Min = -1.0e-3
        y1Max = +2.0e-3
        useLogScale1 = False
        y1Label = ''

    elif plotWhat == 'thermo':

        # Thermodynamics

        useTwoAxes = True

        # Fields to plot
        fields0 = np.array( [ 'AF_T' ] )
        yScale0 = np.array( [ 1.0 ], np.float64 )
        label0 = [ r'$T$' ]
        useCustomY0Limits = False
        y0Min = 6.0e9
        y0Max = 8.0e9
        useLogScale0 = True
        y0Label = r'$\mathrm{K}$'

        fields1 = np.array( [ 'AF_S' ] )
        yScale1 = np.array( [ 1.0 ], np.float64 )
        label1 = [ r'$s$' ]
        useCustomY1Limits = False
        y1Min = 2.0
        y1Max = 14.0
        useLogScale1 = False
        y1Label = r'$k_{\mathrm{B}}/\mathrm{baryon}$'

    elif plotWhat == 'ye':

        # Plot electron fraction only

        useTwoAxes = False

        fields0 = np.array( [ 'AF_Ye' ] )
        yScale0 = np.array( [ 1.0 ], np.float64 )
        label0 = [ r'$Y_{e}$' ]
        useCustomY0Limits = False
        y0Min = -1.0e23
        y0Max = 1.0e1
        useLogScale0 = False
        y0Label = ''

        # Dummy data

        fields1 = np.array( [ 'PF_V1' ] )
        yScale1 = np.array( [ 2.99792458e5 ], np.float64 )
        label1 = [ r'$v^{1}/c$' ]
        useCustomY1Limits = False
        y1Min = -1.0e-3
        y1Max = +2.0e-3
        useLogScale1 = False
        y1Label = ''

    elif plotWhat == 'metric':

        # Metric

        useTwoAxes = True

        # Fields to plot
        fields0 = np.array( [ 'GF_Psi', 'GF_Alpha' ] )
        yScale0 = np.array( [ 1.0, 1.0 ], np.float64 )
        label0 = [ r'$\Psi$', r'$\alpha$' ]
        useCustomY0Limits = False
        y0Min = -1.0e23
        y0Max = 1.0e1
        useLogScale0 = False
        y0Label = ''

        fields1 = np.array( [ 'GF_Beta_1' ] )
        yScale1 = np.array( [ 2.99792458e5 ], np.float64 )
        label1 = [ r'$\beta^{1}/c$' ]
        useCustomY1Limits = False
        y1Min = -1.0e-3
        y1Max = +2.0e-3
        useLogScale1 = False
        y1Label = ''

    else:

        useTwoAxes = False

        fields0 = np.array( [ plotWhat ] )
        yScale0 = np.array( [ 1.0 ], np.float64 )
        label0 = [ r'$T$' ]
        useCustomY0Limits = False
        y0Min = 0.9e10
        y0Max = 1.1e11
        useLogScale0 = True
        y0Label = ''

        # DUMMY
        fields1 = np.array( [ 'GF_Beta_1' ] )
        yScale1 = np.array( [ 2.99792458e5 ], np.float64 )
        label1 = [ r'$\beta^{1}/c$' ]
        useCustomY1Limits = False
        y1Min = -1.0e-3
        y1Max = +2.0e-3
        useLogScale1 = False
        y1Label = ''

    return fields0, yScale0, label0, useCustomY0Limits, \
           y0Min, y0Max, useLogScale0, y0Label, \
           fields1, yScale1, label1, useCustomY1Limits, \
           y1Min, y1Max, useLogScale1, y1Label, useTwoAxes

fields0, yScale0, label0, useCustomY0Limits, \
y0Min, y0Max, useLogScale0, y0Label, \
fields1, yScale1, label1, useCustomY1Limits, \
y1Min, y1Max, useLogScale1, y1Label, useTwoAxes = setFields( 'hydro' )

# Unit system of the data
usePhysicalUnits = True

# Coordinate system (currently supports 'cartesian' and 'spherical' )
coordinateSystem = 'spherical'

# First and last snapshots and number of snapshots to include in movie
SSi = -1 # -1 -> SSi = 0
SSf = -1 # -1 -> SSf = PlotFileArray.shape[0] - 1
nSS = 100 # -1 -> nSS = PlotFileArray.shape[0]

# Max level of refinement to include
maxLevel = -1

plotMesh = False

verbose = True

movieRunTime = 10.0 # seconds

#### ====== End of User Input =======

dataFileDirectory0 = []
for i in range( fields0.shape[0] ):
    dataFileDirectory0.append \
      ( '.{:s}_{:s}_MovieData1D/'.format( ID, fields0[i] ) )
dataFileDirectory0 = np.array( dataFileDirectory0 )

if useTwoAxes:
    dataFileDirectory1 = []
    for i in range( fields1.shape[0] ):
        dataFileDirectory1.append \
          ( '.{:s}_{:s}_MovieData1D/'.format( ID, fields1[i] ) )
    dataFileDirectory1 = np.array( dataFileDirectory1 )

tmpField = ''
for i in range( fields0.shape[0] ):
    tmpField += '_' + fields0[i]
if useTwoAxes :
    for i in range( fields1.shape[0] ):
        tmpField += '_' + fields1[i]

movieName = 'mov.{:s}{:s}.mp4'.format( ID, tmpField )

fig = plt.figure( figsize = (10,6) )
ax0  = fig.add_subplot( 111 )
#ax0.grid( axis = 'x', which = 'both' )

if useCustomY0Limits: ax0.set_ylim( y0Min, y0Max )
if useTwoAxes:
    ax1  = ax0.twinx()
    if useCustomY1Limits: ax1.set_ylim( y1Min, y1Max )

# Append "/" if not present
if plotFileDirectory0[-1] != '/' : plotFileDirectory0 += '/'

if useTwoAxes:
    # Append "/" if not present
    if plotFileDirectory1[-1] != '/' : plotFileDirectory1 += '/'

PlotFileArray0 \
  = GetFileArray \
      ( plotFileDirectory0, plotFileBaseName, \
        SSi = SSi, SSf = SSf, nSS = nSS )

if useTwoAxes:
    PlotFileArray1 \
      = GetFileArray \
          ( plotFileDirectory1, plotFileBaseName, \
            SSi = SSi, SSf = SSf, nSS = nSS )

if nSS < 0: nSS = PlotFileArray0.shape[0]

y0Min = +np.inf
y0Max = -np.inf
for i in range( fields0.shape[0] ):

    MakeDataFile( fields0[i], plotFileDirectory0, dataFileDirectory0[i], \
                  plotFileBaseName, \
                  SSi = SSi, SSf = SSf, nSS = nSS, \
                  maxLevel = maxLevel, verbose = verbose, \
                  forceChoice = False, OW = True )

    if not useCustomY0Limits:
        if verbose: print( 'Getting ylims for {:}...'.format( fields0[i] ) )
        for j in range( nSS ):
            dataFile = dataFileDirectory0[i] + PlotFileArray0[j] + '.dat'
            Data = np.loadtxt( dataFile ) / yScale0[i]
            y0Min = min( y0Min, Data.min() )
            y0Max = max( y0Max, Data.max() )

if not useCustomY0Limits: ax0.set_ylim( y0Min, y0Max )

if useTwoAxes:

    y1Min = +np.inf
    y1Max = -np.inf
    for i in range( fields1.shape[0] ):

        MakeDataFile( fields1[i], plotFileDirectory1, dataFileDirectory1[i], \
                      plotFileBaseName, \
                      SSi = SSi, SSf = SSf, nSS = nSS, \
                      maxLevel = maxLevel, verbose = verbose, \
                      forceChoice = False, OW = True )

        if not useCustomY1Limits:
            if verbose: print( 'Getting ylims for {:}...'.format( fields1[i] ) )
            for j in range( nSS ):
                dataFile = dataFileDirectory1[i] + PlotFileArray1[j] + '.dat'
                Data = np.loadtxt( dataFile ) / yScale1[i]
                y1Min = min( y1Min, Data.min() )
                y1Max = max( y1Max, Data.max() )

    if not useCustomY1Limits: ax1.set_ylim( y1Min, y1Max )

def f(t,i,F=0):

    if t != 0 and i == 0 :
        if t+1 % 10 == 0 :
          print( '{:d}/{:d}'.format( t, nSS ) )

    if F == 0:
        dataFile = dataFileDirectory0[i] + PlotFileArray0[t] + '.dat'
    else:
        dataFile = dataFileDirectory1[i] + PlotFileArray1[t] + '.dat'

    dataShape, dataUnits, time, X1_C, X2_C, X3_C, dX1, dX2, dX3 \
      = ReadHeader( dataFile )

    data = np.loadtxt( dataFile ).reshape( np.int64( dataShape ) )

    if F == 0:
        data /= yScale0[i]
    else:
        data /= yScale1[i]

    return data, dataUnits, X1_C, dX1, time

# this is just to get the mesh info
data0, dataUnits, X1_C, dX1, time = f(0,0)

xL = np.array( [ X1_C[0 ]-0.5*dX1[0 ] ], np.float64 )
xU = np.array( [ X1_C[-1]+0.5*dX1[-1] ], np.float64 )

timeUnits = ''
lengthUnits = ''
if usePhysicalUnits:
    timeUnits   = 'ms'
    lengthUnits = 'km'

ax0.set_xlim( xL[0], xU[0] )

if useTwoAxes:
   ax1.set_xlim( xL[0], xU[0] )

ax0.set_xlabel( 'X1' + ' ' + lengthUnits )

time_text \
  = plt.text( 0.3, 0.9, '', fontsize = 15, transform = ax0.transAxes )

if useLogScale0 : ax0.set_yscale( 'symlog', linthresh = 1.0e20 )
if useTwoAxes:
    if useLogScale1 : ax1.set_yscale( 'log' )

line0 = np.empty( fields0.shape[0], object )
c = [ 'r-', 'rx' ]
for i in range( fields0.shape[0] ):
    line0[i], = ax0.plot( [], [], c[i], label = label0[i] )
if not useTwoAxes:
    ax0.tick_params( axis = 'y', colors = 'k', labelsize = 15 )
else:
    ax0.tick_params( axis = 'y', colors = 'r', labelsize = 15 )
    ax0.set_ylabel( y0Label, c = 'r' )

if useTwoAxes :

    line1 = np.empty( fields1.shape[0], object )
    c = [ 'b-', 'bx' ]
    for i in range( fields1.shape[0] ):
        line1[i], = ax1.plot( [],[], c[i], label = label1[i] )
    ax1.tick_params( axis = 'y', colors = 'b', labelsize = 15 )
    ax1.set_ylabel( y1Label, c = 'b' )

if plotMesh: mesh, = ax0.plot( [], [], 'm.', label = 'mesh' )
ax0.tick_params( axis = 'x', colors = 'k', labelsize = 15 )
ax0.legend( loc = 2, prop = { 'size' : 15 } )
if useTwoAxes:
    ax1.legend( loc = 1, prop = { 'size' : 15 } )

def initializeFrame():

    ret = []

    for i in range( fields0.shape[0] ):
        line0[i].set_data([],[])
        ret.append( line0[i] )

    if useTwoAxes:
      for i in range( fields1.shape[0] ):
          line1[i].set_data([],[])
          ret.append( line1[i] )

    time_text.set_text('')
    ret.append( time_text )

    if plotMesh:
        mesh.set_data([],[])
        ret.append( mesh )

    return ( ret )

def updateFrame( t ):

    ret = []

    for i in range( fields0.shape[0] ):
        data0, dataUnits0, X1_C, dX1, time = f(t,i)
        line0[i].set_data( X1_C, data0.flatten() )
        ret.append( line0[i] )

    if useTwoAxes:
        for i in range( fields1.shape[0] ):
            data1, dataUnits1, X1_C, dX1, time = f(t,i,F=1)
            line1[i].set_data( X1_C, data1.flatten() )
            ret.append( line1[i] )

    if plotMesh:
        mesh.set_data( X1_C - 0.5 * dX1, np.ones( X1_C.shape[0] ) )
        ret.append( mesh )

    time_text.set_text( r'$t-t_b$ = {:.3e} {:}'.format( time-tb, timeUnits ) )
    ret.append( time_text )

    return ( ret )

anim = animation.FuncAnimation( fig, updateFrame, \
                                init_func = initializeFrame, \
                                frames = nSS, \
                                blit = True )

fps = max( 1, nSS // movieRunTime )

print( '\nSaving movie {:}...'.format( movieName ) )
anim.save( movieName, fps = fps, dpi = 300 )
print( 'Done!' )

import os
os.system( 'rm -rf __pycache__ ' )
