#!/usr/bin/env python3

import yt
import numpy as np
import matplotlib.pyplot as plt
import os

from UtilitiesModule import ChoosePlotFile, Overwrite, GetData, GetFileArray

yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

def MakeLineOutPlot( plotFileDirectory, plotFileBaseName, entropyThreshold ):

    Data, DataUnit, r, theta, phi, dr, dtheta, dphi, \
    xL, xU, nX, Time \
      = GetData( plotFileDirectory, plotFileBaseName, \
                 'PolytropicConstant', 'spherical', True, \
                 argv = [ 'a' ], \
                 ReturnMesh = True, ReturnTime = True )

    plt.semilogy( r, Data, 'k-' )
    plt.axhline( entropyThreshold, label = 'Shock radius cut-off' )

    plt.legend()
    plt.savefig( '{:}_EntropyThresholdCheck.png'.format( self.ID ), \
                 dpi = 300 )
    #plt.show()
    plt.close()

    return


def MakeDataFile \
      ( plotFileDirectory, plotFileBaseName, dataFileName, forceChoice, OW ):

    OW = Overwrite( dataFileName, ForceChoice = forceChoice, OW = OW )

    if not OW: return

    plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )

    # Just to get number of elements
    plotFile = plotFileDirectory + plotFileArray[0]
    dum, dumm, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile, 'PolytropicConstant' )

    if( nX[1] > 1 ):
        Data = np.empty( (plotFileArray.shape[0],nX[0],nX[1]), np.float64 )
    else:
        Data = np.empty( (plotFileArray.shape[0],nX[0]), np.float64 )

    Time = np.empty( (plotFileArray.shape[0]), np.float64 )

    for i in range( plotFileArray.shape[0] ):

        plotFile = plotFileDirectory + plotFileArray[i]

        Time[i], Data[i], dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
          = GetData( plotFile, 'PolytropicConstant' )

    # Save multi-D array with np.savetxt. Taken from:
    # https://stackoverflow.com/questions/3685265/
    # how-to-write-a-multidimensional-array-to-a-text-file
    with open( dataFileName, 'w' ) as FileOut:

        FileOut.write( '# Array shape: {0}\n'.format( Data.shape ) )

        # Iterating through an n-dimensional array produces slices along
        # the last axis. This is equivalent to Data[i] in this case
        i = 0
        for TimeSlice in Data:

            i += 1

            np.savetxt( FileOut, TimeSlice )
            FileOut.write( '# New slice\n' )

    return


def ComputeShockRadius \
       ( shockRadiusFileName, forceChoice, OW, dataFileName, \
         plotFileDirectory, plotFileBaseName, entropyThreshold ):

    OW = Overwrite( shockRadiusFileName, forceChoice, OW )

    if not OW: return

    MakeLineOutPlot( plotFileDirectory, plotFileBaseName, entropyThreshold )

    plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )

    # Just to get number of elements
    plotFile = plotFileDirectory + plotFileArray[0]
    dummy, dummy2, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile, 'PolytropicConstant' )
    dX1 = dX[0]
    dX2 = dX[1]

    OW = Overwrite( dataFileName )

    MakeDataFile \
      ( plotFileDirectory, plotFileBaseName, dataFileName, forceChoice, OW )

    Data = np.loadtxt( dataFileName ).reshape( \
             (plotFileArray.shape[0],nX[0],nX[1]) )
    Time = np.loadtxt( TimeFileName )

    DataUnit = 'erg/cm**3/(g/cm**3)**(4/3)'

    Volume = np.zeros( PlotFileArray.shape[0], np.float64 )
    Min    = np.empty( PlotFileArray.shape[0], np.float64 )
    Max    = np.empty( PlotFileArray.shape[0], np.float64 )

    if self.Verbose:
        print( '\n  Computing shock radius...' )
        print(     '-------------------------' )

    for iT in range( Data.shape[0] ):

        if self.Verbose and ( (iT+1) % 10 == 0 ):
            print( '    File {:}/{:}'.format( iT+1, Data.shape[0] ) )

        # Sum volumes of elements containing shocked fluid

        ind = np.where( Data[iT] < self.EntropyThreshold )
        Min[iT] = X1[ind[0].min()]

        ind = np.where( Data[iT] > self.EntropyThreshold )
        Max[iT] = X1[ind[0].max()]

        for iX1 in range( Data.shape[1] ):

            for iX2 in range( Data.shape[2] ):

                if( Data[iT,iX1,iX2] > self.EntropyThreshold ):

                    Volume[iT] \
                      += 2.0 * np.pi \
                           * 1.0 / 3.0 * ( ( X1[iX1] + dX1 )**3 \
                                               - X1[iX1]**3 ) \
                           * ( np.cos( X2[iX2] ) - np.cos( X2[iX2] + dX2 ) )

        Rs = ( Volume[iT] / ( 4.0 / 3.0 * np.pi ) )**( 1.0 / 3.0 )
        if Min[iT] >= Rs: Min[iT] = Rs
        if Max[iT] <= Rs: Max[iT] = Rs

    # 4/3 * pi * Rs^3 = Volume
    # ==> Rs = ( Volume / ( 4/3 * pi ) )^( 1 / 3 )
    ShockRadius = ( Volume / ( 4.0 / 3.0 * np.pi ) )**( 1.0 / 3.0 )

    np.savetxt( self.ShockRadiusVsTimeFileName, \
                np.vstack( ( Time, ShockRadius, Min, Max ) ) )

    os.system( 'rm -f {:}'.format( self.TimeFileName  ) )
    os.system( 'rm -f {:}'.format( self.DataFileName  ) )

    return


if __name__ == "__main__":

    RootDirectory = '/lump/data/accretionShockStudy/'
    EntropyThreshold = 1.0e15

    SaveFigAs                 \
      = 'fig.ShockRadiusVsTime_CompareGamma_{:}.png'.format \
        ( 'PolytropicConstant' )

    plt.suptitle( 'GR1D_M2.0_Mdot0.3_Rs150' )

    ID = 'GR1D_M2.0_Mdot0.3_Rs150'
    SR = ShockRadius( RootDirectory, ID, EntropyThreshold = EntropyThreshold, \
                      PlotFileDirectorySuffix = '', \
                      ForceChoice = False, OW = True, ModPlotFile = 1, \
                      Verbose = True )
    SR.ComputeShockRadius()
    t133, r133, Min133, Max133 = np.loadtxt( SR.ShockRadiusVsTimeFileName )
    r133_0 = r133[0]
    plt.plot( t133, ( r133 - r133_0 ) / r133_0, 'r.', label = r'$\Gamma=4/3$' )

    ID = 'GR1D_M2.0_Mdot0.3_Rs150_Gm1.13_nX2560'
    SR = ShockRadius( RootDirectory, ID, EntropyThreshold = EntropyThreshold, \
                      PlotFileDirectorySuffix = '', \
                      ForceChoice = False, OW = True, ModPlotFile = 1, \
                      Verbose = True )
    SR.ComputeShockRadius()
    t113, r113, Min113, Max113 = np.loadtxt( SR.ShockRadiusVsTimeFileName )
    r113_0 = r113[0]

    plt.plot( t113, ( r113 - r113_0 ) / r113_0, 'b.', label = r'$\Gamma=1.13$' )

    plt.xlim( min( t133[0], t113[0] ), max( t133[-1], t113[-1] ))

    plt.xlabel( 'Time [ms]' )
    plt.ylabel( r'$\left(R_{s}-R_{s,0}\right)/R_{s,0}$', labelpad = -0.05 )
    plt.legend()
    #plt.show()
    plt.savefig( SaveFigAs, dpi = 300 )
    plt.close()

    import os
    os.system( 'rm -rf __pycache__ ' )
