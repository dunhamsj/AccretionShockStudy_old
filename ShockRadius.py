#!/usr/bin/env python3

import yt
import numpy as np
import matplotlib.pyplot as plt
import os

from UtilitiesModule import ChoosePlotFile, Overwrite, GetData, GetFileArray

yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

class ShockRadius:

    def __init__( self, RootDirectory, ID, EntropyThreshold = 1.0e15, \
                  PlotFileDirectorySuffix = '', \
                  ForceChoice = False, OW = True, ModPlotFile = 1, \
                  Verbose = False ):


        self.RootDirectory           = RootDirectory
        self.ID                      = ID
        self.EntropyThreshold        = EntropyThreshold
        self.PlotFileDirectorySuffix = PlotFileDirectorySuffix
        self.ForceChoice             = ForceChoice
        self.OW                      = OW
        self.ModPlotFile             = ModPlotFile
        self.Verbose                 = Verbose

        if self.Verbose:

            print( '\n  Creating instance of ShockRadius class...\n' )
            print( '    Variables:' )
            print( '    ----------' )
            print( '{:>29}: {:}'.format( 'RootDirectory', self.RootDirectory ) )
            print( '{:>29}: {:}'.format( 'ID', self.ID ) )
            print( '{:>29}: {:.3e} cgs'.format \
              ( 'EntropyThreshold', self.EntropyThreshold ) )
            print( '{:>29}: {:}'.format \
              ( 'PlotFileDirectorySuffix', self.PlotFileDirectorySuffix ) )
            print( '{:>29}: {:}'.format \
              ( 'ForceChoice', self.ForceChoice ) )
            print( '{:>29}: {:}'.format \
              ( 'OW', self.OW ) )
            print( '{:>29}: {:}'.format \
              ( 'ModPlotFile', self.ModPlotFile ) )
            print('')

        self.DataFileName              \
          = '.{:}_ShockRadiusVsTimeData_{:}.dat'.format \
            ( ID, 'PolytropicConstant' )
        self.TimeFileName              \
          = '.{:}_ShockRadiusVsTimeTime.dat'.format \
            ( ID )
        self.ShockRadiusVsTimeFileName \
          = '.{:}_ShockRadiusVsTime.dat'.format \
            ( ID )
        self.SaveFigAs                 \
          = 'fig.ShockRadiusVsTime_{:}_{:}.png'.format \
            ( ID, 'PolytropicConstant' )

        self.PlotFileDirectory \
          = self.RootDirectory \
              + self.ID + '{:}/'.format( PlotFileDirectorySuffix )

        self.PlotFileBaseName = self.ID + '.plt'

        self.PlotFileArray \
          = GetFileArray( self.PlotFileDirectory, self.PlotFileBaseName )

        self.PlotFileArray = self.PlotFileArray[::self.ModPlotFile]

        File \
          = ChoosePlotFile( self.PlotFileDirectory, self.PlotFileBaseName )

        self.ds = yt.load( '{:}'.format( self.PlotFileDirectory + File ) )

        self.nX = self.ds.domain_dimensions
        self.xL = self.ds.domain_left_edge.to_ndarray()
        self.xH = self.ds.domain_right_edge.to_ndarray()

        self.dX = ( self.xH - self.xL ) / np.float64( self.nX )

        self.X1 = np.linspace( self.xL[0] + 0.5 * self.dX[0], \
                               self.xH[0] - 0.5 * self.dX[0], \
                               self.nX[0] )

        self.X2 = np.linspace( self.xL[1] + 0.5 * self.dX[1], \
                               self.xH[1] - 0.5 * self.dX[1], \
                               self.nX[1] )

        return


    def MakeLineOutPlot( self ):

        Data, DataUnit, r, theta, phi, dr, dtheta, dphi, \
        xL, xU, nX, Time \
          = GetData( self.PlotFileDirectory, self.PlotFileBaseName, \
                     'PolytropicConstant', 'spherical', True, \
                     argv = [ 'a' ], \
                     ReturnMesh = True, ReturnTime = True )

        fig2 = plt.figure()

        plt.semilogy( r, Data, 'k-' )
        plt.axhline( self.EntropyThreshold, label = 'Shock radius cut-off' )

        plt.legend()
        plt.savefig( '{:}_EntropyThresholdCheck.png'.format( self.ID ), \
                     dpi = 300 )
        #plt.show()
        plt.close()

        del fig2

        input( '\nPress enter to proceed...' )

        return


    def MakeDataFile( self ):

        OW = Overwrite( self.DataFileName, \
                        ForceChoice = self.ForceChoice, OW = self.OW )

        if not OW: return

        PlotFileArray     = self.PlotFileArray
        nX                = self.nX
        PlotFileDirectory = self.PlotFileDirectory
        PlotFileBaseName  = self.PlotFileBaseName
        DataFileName      = self.DataFileName
        TimeFileName      = self.TimeFileName

        # Put all time-slices into one array to use for movie making

        if nX[1] > 1:
            Data = np.empty( (PlotFileArray.shape[0],nX[0],nX[1]), np.float64 )
        else:
            Data = np.empty( (PlotFileArray.shape[0],nX[0]), np.float64 )

        Time = np.empty( (PlotFileArray.shape[0]), np.float64 )

        if self.Verbose:
            print( '\n    Making data array for shock radius vs. time...' )
            print(   '    ----------------------------------------------' )

        for i in range( PlotFileArray.shape[0] ):

            if self.Verbose and ( (i+1) % 10 == 0 ):
                print( '      File {:}/{:}'.format \
                  ( i+1, PlotFileArray.shape[0] ) )

            Data[i], DataUnit, r, theta, phi, dr, dtheta, dphi, \
            xL, xU, nX, Time[i] \
              = GetData( PlotFileDirectory, PlotFileBaseName, \
                         'PolytropicConstant', 'spherical', True, \
                         argv = [ 'a', PlotFileArray[i] ], \
                         ReturnMesh = True, ReturnTime = True )

        np.savetxt( TimeFileName, Time )

        if self.Verbose:
            print( '\n  Making data file for shock radius vs. time...' )
            print(   '---------------------------------------------' )

        # Save multi-D array with np.savetxt. Taken from:
        # https://stackoverflow.com/questions/3685265/
        # how-to-write-a-multidimensional-array-to-a-text-file
        with open( DataFileName, 'w' ) as FileOut:

            FileOut.write( '# Array shape: {0}\n'.format( Data.shape ) )

            # Iterating through an n-dimensional array produces slices along
            # the last axis. This is equivalent to Data[i] in this case
            i = 0
            for TimeSlice in Data:

                i += 1

                if self.Verbose and ( i % 10 == 0 ):
                    print( '    File {:}/{:}'.format \
                      ( i, PlotFileArray.shape[0] ) )

                np.savetxt( FileOut, TimeSlice )
                FileOut.write( '# New slice\n' )

        return


    def ComputeShockRadius( self ):

        OW = Overwrite( self.ShockRadiusVsTimeFileName, \
                        ForceChoice = self.ForceChoice, OW = self.OW )

        if not OW: return

        self.MakeLineOutPlot()

        if self.Verbose:
            print( '\n    Computing average shock radius...' )

        PlotFileArray = self.PlotFileArray
        nX            = self.nX
        DataFileName  = self.DataFileName
        TimeFileName  = self.TimeFileName
        X1            = self.X1
        X2            = self.X2
        dX            = self.dX

        dX1 = dX[0]
        dX2 = dX[1]

        OW = Overwrite( self.DataFileName, \
                        ForceChoice = self.ForceChoice, OW = self.OW )

        self.MakeDataFile()

        Data = np.loadtxt( DataFileName ).reshape( \
                 (PlotFileArray.shape[0],nX[0],nX[1]) )
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

    def PlotShockRadiusVsTime( self ):

        t, r, Min, Max = np.loadtxt( self.ShockRadiusVsTimeFileName )

        r0 = r[0]

        plt.plot( t, ( r - r0 ) / r0, 'k-', label = 'Average' )
        #plt.plot( t, Max/r0, 'b-', label = 'Max' )
        #plt.plot( t, Min/r0, 'r-', label = 'Min' )
        plt.xlim( t[0], t[-1] )
        plt.xlabel( 'Time [ms]' )
        plt.ylabel( r'$\left(R_{s}-R_{s,0}\right)/R_{s,0}$' )
        plt.legend()
        #plt.show()
        plt.savefig( self.SaveFigAs, dpi = 300 )
        plt.close()


if __name__ == "__main__":

    #RootDirectory = '/Users/dunhamsj/Work/Data/AccretionShockParameterStudy/'
    RootDirectory = '/lump/data/AccretionShockStudy/'
    EntropyThreshold = 1.0e15

    ID = 'GR1D_M2.4_Mdot0.3_Rs165'

    SR = ShockRadius( RootDirectory, ID, EntropyThreshold = EntropyThreshold, \
                      PlotFileDirectorySuffix = '', \
                      ForceChoice = False, OW = True, ModPlotFile = 10, \
                      Verbose = True )

    SR.ComputeShockRadius()
    SR.PlotShockRadiusVsTime()

    import os
    os.system( 'rm -rf __pycache__ ' )
