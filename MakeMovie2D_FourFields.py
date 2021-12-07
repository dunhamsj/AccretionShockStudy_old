#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import LogNorm, SymLogNorm
import subprocess
import os

from MakeDataFile import MakeDataFile
from UtilitiesModule import GetNorm


def MakeMovie2D( ID, Field, MovieName, Title, Labels ):

    """
    Generate a movie from data files created in MakeDataFile.py.
    """

    # --- Get user's HOME directory ---

    HOME = subprocess.check_output( ["echo $HOME"], shell = True)
    HOME = HOME[:-1].decode( "utf-8" ) + '/'

    nPanels = ID.shape[0]

    cmap           = np.empty( nPanels, object )
    UseLogScale    = np.empty( nPanels, bool )
    UseCustomTicks = np.empty( nPanels, bool )

    ############# User input #############

    ccmap = 'RdBu'
    ULS   = True
    UCT   = False

    cmap          [0] = ccmap
    UseLogScale   [0] = ULS
    UseCustomTicks[0] = UCT

    if( nPanels > 1 ):

        cmap          [1] = ccmap
        UseLogScale   [1] = ULS
        UseCustomTicks[1] = UCT

    if( nPanels > 2 ):

        cmap          [2] = ccmap
        UseLogScale   [2] = ULS
        UseCustomTicks[2] = UCT

    if( nPanels > 3 ):

        cmap          [3] = ccmap
        UseLogScale   [3] = ULS
        UseCustomTicks[3] = UCT

    MovieRunTime = 10.0 # seconds

    ############# End of user input #############

    Half  = 1.0 / 2.0
    Third = 1.0 / 3.0

    DataDirectory    = np.empty( nPanels, object )
    PlotFileBaseName = np.empty( nPanels, object )
    DataFileName     = np.empty( nPanels, object )
    Data             = np.empty( nPanels, object )
    X1               = np.empty( nPanels, object )
    X2               = np.empty( nPanels, object )
    dr               = np.empty( nPanels, np.float64 )
    dt               = np.empty( nPanels, np.float64 )
    xL               = np.empty( (nPanels,3), np.float64 )
    xU               = np.empty( (nPanels,3), np.float64 )
    nX               = np.empty( (nPanels,3), np.int64 )
    FieldUnit        = np.empty( nPanels, object )
    DataShape        = np.empty( (nPanels,3), np.int64 )
    vmin             = np.empty( nPanels, np.float64 )
    vmax             = np.empty( nPanels, np.float64 )
    nTicks           = np.empty( nPanels, np.int64 )
    ticks            = np.empty( nPanels, object )
    ticklabels       = np.empty( nPanels, object )
    Norm             = np.empty( nPanels, object )
    r                = np.empty( nPanels, object )
    theta            = np.empty( nPanels, object )
    im               = np.empty( nPanels, object )
    f                = np.empty( nPanels, object )
    cbar             = np.empty( nPanels, object )
    cbaxes           = np.empty( nPanels, object )
    IDi              = np.empty( nPanels, object )

    fig = plt.figure( figsize = (16,9) )
    fig.suptitle( Title )
    ax  = fig.add_subplot( 111, polar = True )

    nFrames = +np.inf

    tMin = +np.inf
    tMax = +np.inf

    for i in range( nPanels ):

        DataDirectory[i] \
          = HOME + 'Research/Data/AccretionShockParameterStudy/{:}/'.format \
                         ( ID[i] )

        PlotFileBaseName[i] = ID[i] + '.plt'
        IDi[i]              = '{:}_{:}'.format( ID[i], Field[i] )
        DataFileName    [i] = '.{:}.dat'.format( IDi[i] )

        xLi, xUi, nXi, FileArray \
          = MakeDataFile( Field[i], DataDirectory[i], \
                          DataFileName[i], PlotFileBaseName[i] )

        xL[i] = xLi.to_ndarray()
        xU[i] = xUi.to_ndarray()
        nX[i] = nXi
        dr[i] = ( xU[i,0] - xL[i,0] ) / np.float64( nX[i,0] )
        dt[i] = ( xU[i,1] - xL[i,1] ) / np.float64( nX[i,1] )

        f            = open( DataFileName[i] )
        header       = f.readline()[16:-2]
        FieldUnit[i] = f.readline()[9:-1]
        DataShape[i] \
          = tuple( [ np.int64( dim ) for dim in header.split( ',' ) ] )

        t = list( [ np.float64( t ) for t in f.readline()[12:-2].split(' ') ] )
        t = np.array( t )

        nFrames = min( nFrames, t.shape[0] )

        tMin = min( tMin, t[0] )
        tMax = min( tMax, t[nFrames-1] )

        print( 'Reading in data...' )

        d = np.loadtxt( DataFileName[i], \
                        dtype = np.float64 ).reshape( DataShape[i] )

        X1[i] \
          = np.linspace( xL[i,0] + Half * dr[i], \
                         xU[i,0] - Half * dr[i], nX[i,0] )

        nT = np.int64( Half * nX[i,1] )

        if   i == 0:

            Data[i] = d[:,:,:nT]

            thetaL = xL[i,1] + Half * dt[i]

        elif i == 1:

            Data[i] = d[:,:,nT:]

            thetaL = xU[i,1] - Half * ( nX[i,1] - 1 ) * dt[i]

        elif i == 2:

            Data[i] = d[:,:,nT:]
            Data[i] = Data[i][:,:,::-1]

            thetaL = xU[i,1] + Half * dt[i]

        elif i == 3:

            Data[i] = d[:,:,:nT]
            Data[i] = Data[i][:,:,::-1]

            thetaL = 2.0 * np.pi - Half * ( nX[i,1] - 1 ) * dt[i]

        thetaU = thetaL + ( nT - 1 ) * dt[i]

        X2[i]  = np.linspace( thetaL, thetaU, nT )
        theta[i], r[i] = np.meshgrid( X2[i], X1[i] )

        if UseLogScale[i]:

            ind = np.where( np.abs( Data[i] ) < 1.0e-16 )
            nZeros = ind[0].shape[0]
            for j in range( nZeros ):
                Data[i][ind[0][j],ind[1][j],ind[2][j]] = 1.0e-17

        vmin[i] = min( +np.inf, Data[i].min() )
        vmax[i] = max( -np.inf, Data[i].max() )

        if UseCustomTicks[i]:

            vMi = 3.0e15
            vMa = 4.6e15

            if i == 0:

              vmin[i] = vMi
              vmax[i] = vMa

            if i == 1:

              vmin[i] = vMi
              vmax[i] = vMa

            if i == 2:

              vmin[i] = vMi
              vmax[i] = vMa

            if i == 3:

              vmin[i] = vMi
              vmax[i] = vMa

            # Custom ticks with SymLog not working yet
            if( UseLogScale[i] and vmin[i] < 0.0 ): UseCustomTicks[i] = False

        if UseLogScale[i]:

            if  ( np.abs( vmax[i] ) / np.abs( vmin[i] ) > 1.0e10 ):
                vmin[i] = vmax[i] / 1.0e10
            elif( np.abs( vmin[i] ) / np.abs( vmax[i] ) > 1.0e10 ):
                vmax[i] = vmin[i] / 1.0e10

    if( np.all( Field == Field[0] ) ):

        VMIN = vmin.min()
        VMAX = vmax.max()

        vmin[:] = VMIN
        vmax[:] = VMAX

    for i in range( nPanels ):

        if( UseCustomTicks[i] ):

            nTicks[i] = 5

            if( UseLogScale[i] ):

                ticks[i] \
                  = np.logspace( np.log10( vmin[i] ), \
                                 np.log10( vmax[i] ), nTicks[i] )

            else:

                ticks[i] = np.linspace( vmin[i], vmax[i], nTicks[i] )

            ticklabels[i] = []
            for tick in ticks[i]:
                ticklabels[i].append( '{:.3e}'.format( tick ) )

    def f(t,i):
        return Data[i][t]

    for i in range( nPanels ):

        Norm[i] = GetNorm( UseLogScale[i], Data[i], \
                           vmin = vmin[i], vmax = vmax[i] )

        im[i] = ax.pcolormesh( theta[i], r[i], f(0,i)[:,:], \
                               cmap = cmap[i], \
                               norm = Norm[i], \
                               shading = 'nearest' )

    # Limits on coordinate axes

    ax.set_thetamin( 0.0 )
    ax.set_thetamax( 360.0 )
    ax.set_theta_direction( -1 )

    ax.set_theta_zero_location( 'N' ) # z-axis vertical
    time_text = plt.text( 0.125 * np.pi / 2.0, xU[0,0] * ( 1.0 + 0.1 ), '' )

    cbaxes[0] = fig.add_axes( [ 0.80, 0.53, 0.03, 0.4 ] )
    if( nPanels > 1 ): cbaxes[1] = fig.add_axes( [ 0.80, 0.05, 0.03, 0.4 ] )
    if( nPanels > 2 ): cbaxes[2] = fig.add_axes( [ 0.20, 0.05, 0.03, 0.4 ] )
    if( nPanels > 3 ): cbaxes[3] = fig.add_axes( [ 0.20, 0.53, 0.03, 0.4 ] )

    for i in range( nPanels ):

        cbar[i] = fig.colorbar( im[i], cax = cbaxes[i] )

        if( UseCustomTicks[i] ):

            cbar[i].ax.minorticks_off()
            cbar[i].set_ticks( ticks[i] )
            cbar[i].ax.set_yticklabels( ticklabels[i] )

        cbar[i].set_label( Labels[i] )#+ ' [' + FieldUnit[i] + ']' )

    if( nPanels > 2 ):

        cbaxes[2].yaxis.set_ticks_position( 'left' )
        cbaxes[2].yaxis.set_label_position( 'left' )

    if( nPanels > 3 ):

        cbaxes[3].yaxis.set_ticks_position( 'left' )
        cbaxes[3].yaxis.set_label_position( 'left' )

    Time = np.linspace( tMin, tMax, nFrames )

    def UpdateFrame(t):
        ret = []
        for i in range( nPanels ):
            im[i].set_array( f(t,i)[:,:].flatten() )
            ret.append( im[i] )
        time_text.set_text( 'time = {:d} ms'.format( np.int( Time[t] ) ) )
        ret.append( time_text )
        ret = tuple( ret )
        return ret

    # Call the animator

    print( 'Making movie...' )

    fps = nFrames / MovieRunTime

    anim \
      = animation.FuncAnimation \
          ( fig, UpdateFrame, frames = nFrames, blit = True )

    anim.save( MovieName, fps = fps )
    plt.close()

    os.system( 'rm -rf __pycache__ ' )

Field = np.array( [ 'DivV2', \
                    'DivV2', \
                    'DivV2' ] )

ID = np.array( [ 'GR2D_M1.4_Mdot0.3_Rs180_PA0.000_nX320x064_HLL', \
                 'GR2D_M1.4_Mdot0.3_Rs180_PA0.000_nX320x064_HLLC', \
                 'GR2D_M1.4_Mdot0.3_Rs180_PA0.000_nX320x064_HYBRID' ] )

Title = Field[0]
Labels = np.array( [ ID[0][-9:], ID[1][-9:], ID[2][-9:]  ] )

MovieName \
  = 'mov.ThreeFields_RiemannSolvers_PA0.000_{:}.mp4'.format( Field[0] )

MakeMovie2D( ID, Field, MovieName, Title, Labels )
