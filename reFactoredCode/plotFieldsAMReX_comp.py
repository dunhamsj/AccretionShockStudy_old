#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from UtilitiesModule import GetData, GetNorm, MapCenterToCorners

"""

Default usage, plots last Plotfile in PlotfileDirectory:

  $ python3 plotFieldsAMReX.py

Alernate usage, plot specific file in PlotfileDirectory:

  $ python3 plotFieldsAMReX.py 10

  will plot the *00000010 Plotfile

"""

#### ========== User Input ==========

# Specify name of problem
ProblemName = [ 'NR2D_M1.4_Mdot0.3_Rs180', \
                'NR2D_M1.4_Rpns040_Rs180_Mdot0.3' ]

# Specify title of figure
FigTitle = ProblemName[1]

# Specify directory containing amrex Plotfiles
PlotfileDirectory \
  = [ '/lump/data/accretionShockStudy/{:}/' \
      .format( ProblemName[0] ), \
      '/lump/data/accretionShockStudy/newRuns/newProductionRuns/{:}/' \
      .format( ProblemName[1] ) ]

# Specify plot file base name
PlotfileBaseName = [ ProblemName[0] + '.plt_', ProblemName[1] + '.plt' ]

# Specify field to plot
Field = 'PF_D'

# Specify to plot in log-scale
UseLogScale_X  = False
UseLogScale_Y  = False
UseLogScale_2D = True

# Specify whether or not to use physical units
UsePhysicalUnits = True

# Specify coordinate system (currently supports 'cartesian' and 'spherical')
CoordinateSystem = 'spherical'

# Max level of refinement to plot (-1 plots leaf elements)
MaxLevel = -1

# Write extra info to screen
Verbose = True

# Use custom limts for y-axis (1D) or colorbar (2D)
UseCustomLimits = False
vmin = 0.0
vmax = 2.0

# Save figure (True) or plot figure (False)
SaveFig = False

# Specify colormap (2D only)
cmap = 'viridis'

#### ====== End of User Input =======

polar = True

TimeUnits = 'ms'
X1Units = 'km'
X2Units = 'rad'

Data00, DataUnit0, X1_C0, X2_C0, X3_C0, dX10, dX20, dX30, xL0, xH0, nX0, Time0 \
  = GetData( PlotfileDirectory[0], PlotfileBaseName[0], Field, \
             CoordinateSystem, UsePhysicalUnits, argv = ['a','0'], \
             MaxLevel = MaxLevel, \
             ReturnTime = True, ReturnMesh = True, Verbose = True )
Data0, DataUnit0, X1_C0, X2_C0, X3_C0, dX10, dX20, dX30, xL0, xH0, nX0, Time0 \
  = GetData( PlotfileDirectory[0], PlotfileBaseName[0], Field, \
             CoordinateSystem, UsePhysicalUnits, argv = ['a','0'], \
             MaxLevel = MaxLevel, \
             ReturnTime = True, ReturnMesh = True, Verbose = True )
#Data0 = ( Data00 - Data0 ) / Data0

Data10, DataUnit1, X1_C1, X2_C1, X3_C1, dX11, dX21, dX31, xL1, xH1, nX1, Time1 \
  = GetData( PlotfileDirectory[1], PlotfileBaseName[1], Field, \
             CoordinateSystem, UsePhysicalUnits, argv = ['a','0'], \
             MaxLevel = MaxLevel, \
             ReturnTime = True, ReturnMesh = True, Verbose = True )
Data1, DataUnit1, X1_C1, X2_C1, X3_C1, dX11, dX21, dX31, xL1, xH1, nX1, Time1 \
  = GetData( PlotfileDirectory[1], PlotfileBaseName[1], Field, \
             CoordinateSystem, UsePhysicalUnits, argv = ['a','0'], \
             MaxLevel = MaxLevel, \
             ReturnTime = True, ReturnMesh = True, Verbose = True )
#Data1 = ( Data10 - Data1 ) / Data1

nDims = 1
if nX0[1] > 1: nDims += 1

if nDims == 1:

    dX10  = np.copy( dX10 [:,0,0] )
    X1_C0 = np.copy( X1_C0[:,0,0] )
    Data0 = np.copy( Data0[:,0,0] )

    dX11  = np.copy( dX11 [:,0,0] )
    X1_C1 = np.copy( X1_C1[:,0,0] )
    Data1 = np.copy( Data1[:,0,0] )

    fig, ax = plt.subplots( 1, 1 )

    ax.plot( X1_C0, Data0, 'k.', \
             label = '{:}'.format( ProblemName[0][-13:-8] ) )
    ax.plot( X1_C1, Data1, 'r.', \
             label = '{:}'.format( ProblemName[1][-13:-8] ) )
    ax.legend()
    if UseLogScale_Y: ax.set_yscale( 'log' )
    if UseCustomLimits: ax.set_ylim( vmin, vmax )
    ax.set_xlim( min( xL0[0], xL1[0] ), max( xH0[0], xH1[0] ) )
    ax.set_xlabel \
      ( r'$x^{{1}}\ \left[\mathrm{{{:}}}\right]$'.format( X1Units ), \
        fontsize = 15 )
    ax.set_ylabel( Field + ' ' + r'$\mathrm{{{:}}}$'.format( DataUnit0 ) )
    ax.grid()

elif nDims == 2:

    X1_C0 = np.copy( X1_C0[:,:,0] )
    X2_C0 = np.copy( X2_C0[:,:,0] )
    dX10  = np.copy( dX10 [:,:,0] )
    dX20  = np.copy( dX20 [:,:,0] )
    Data0 = np.copy( Data0[:,:,0] )

    X1_C1 = np.copy( X1_C1[:,:,0] )
    X2_C1 = np.copy( X2_C1[:,:,0] )
    dX11  = np.copy( dX11 [:,:,0] )
    dX21  = np.copy( dX21 [:,:,0] )
    Data1 = np.copy( Data1[:,:,0] )

    '''
    # Line-out plots
    fig, ax = plt.subplots( 1, 1 )

    nX2 = X1_C.shape[1]

    ax.plot( X1_C[:,0      ], Data[:,0     ], \
             label = r'$\theta/\pi={:.3f}$'.format( X2_C[0,0     ] / np.pi ) )
    ax.plot( X1_C[:,nX2//2 ], Data[:,nX2//2], \
              label = r'$\theta/\pi={:.3f}$'.format( X2_C[0,nX2//2] / np.pi ) )
    ax.plot( X1_C[:,-1     ], Data[:,-1    ], \
              label = r'$\theta/\pi={:.3f}$'.format( X2_C[0,-1    ] / np.pi ) )

    if UseLogScale_X: ax.set_xscale( 'log' )
    if UseLogScale_Y: ax.set_yscale( 'log' )
    ax.legend()

    plt.show()
    exit()
    '''

    if not UseCustomLimits:

        vmin = min( Data0.min(), Data1.min() )
        vmax = max( Data0.max(), Data1.max() )

    Norm = GetNorm( UseLogScale_2D, Data0, vmin = vmin, vmax = vmax )

    fig = plt.figure( figsize = (12,9) )
    ax = fig.add_subplot( 111, polar = polar )

    # pcolormesh wants the corners of the elements
    X10c, X20c = MapCenterToCorners( X1_C0, X2_C0, dX10, dX20 )
    X11c, X21c = MapCenterToCorners( X1_C1, X2_C1, dX11, dX21 )

    ax.grid( False )
    ax.set_thetamin( 0.0  )
    ax.set_thetamax( 360.0)
    ax.set_rmin( 0.0  )
    ax.set_rmax( max( xH0[0], xH1[1] ) )
    ax.set_theta_direction( -1 )
    ax.set_theta_zero_location( 'N' ) # z-axis vertical

#    ax.text( 0.6, 0.9, r'$t={:.2e}\ \left[\mathrm{{{:}}}\right]$'.format \
#             ( Time0, TimeUnits ), transform = ax.transAxes )

    # Transpose data for spherical-polar coordinates

    X10c = np.copy( X10c[:,0] )
    X20c = np.copy( X20c[0,:] )

    X11c = np.copy( X11c[:,0] )
    X21c = 2.0 * np.pi - np.copy( X21c[0,:] )

    X10c, X20c = np.meshgrid( X20c, X10c )
    X11c, X21c = np.meshgrid( X21c, X11c )

    im0 = ax.pcolormesh( X10c, X20c, Data0, \
                         cmap = cmap, \
                         norm = Norm, \
                         shading = 'flat' )

    im1 = ax.pcolormesh( X11c, X21c, Data1, \
                         cmap = cmap, \
                         norm = Norm, \
                         shading = 'flat' )

    cbar0 = fig.colorbar( im0 )
    cbar1 = fig.colorbar( im1, location = 'left' )

    cbar0.set_label( ProblemName[0][-13:-8], labelpad = -10 )
    cbar1.set_label( ProblemName[1][-13:-8] )

ax.set_title( r'$\texttt{{{:}}}$'.format( FigTitle ) )

if SaveFig:

    FigName = '/home/kkadoogan/fig.IC_1D_{:}.png'.format( Field )
    plt.savefig( FigName, dpi = 300 )

else:

    plt.show()

plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
