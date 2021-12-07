#!/usr/bin/env python3

import yt
import numpy as np
import matplotlib.pyplot as plt
import os

from UtilitiesModule import ChoosePlotFile, OverwriteFile, GetData, GetFileArray

yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

class ShockRadius:

    def __init__( self, Root, ID, Field, EntropyThreshold = 1.0e15 ):

        print( '\nCreating instance of ShockRadius class...\n' )

        self.Root = Root
        self.ID = ID
        self.Field = Field
        self.EntropyThreshold = EntropyThreshold

        print( '  Variables:' )
        print( '  ----------' )
        print( '    Root:             {:s}'.format( self.Root ) )
        print( '    ID:               {:s}'.format( self.ID ) )
        print( '    Field:            {:s}'.format( self.Field ) )
        print( '    EntropyThreshold: {:.3e}\n'.format \
          ( self.EntropyThreshold ) )

        self.DataFileName              \
          = '.{:}_ShockRadiusVsTimeData_{:}.dat'.format( ID, Field )
        self.TimeFileName              \
          = '.{:}_ShockRadiusVsTimeTime.dat'.format( ID )
        self.ShockRadiusVsTimeFileName \
          = '.{:}_ShockRadiusVsTime.dat'.format( ID )
        self.SaveFigAs                 \
          = 'fig.ShockRadiusVsTime_{:}_{:}.png'.format( ID, Field )

        self.DataDirectory = self.Root + '/'

        self.PlotFileBaseName = self.ID + '.plt'

        self.FileArray \
          = GetFileArray( self.DataDirectory, self.PlotFileBaseName )

        File \
          = ChoosePlotFile( self.DataDirectory, self.PlotFileBaseName )

        self.ds = yt.load( '{:}'.format( self.DataDirectory + File ) )

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

        print( '\nCalling ShockRadius.MakeLineOutPlot...\n' )

        ds = self.ds

        X1 = self.X1

        oray = ds.ortho_ray( axis = 0, coords = (0,0) )

        PolytropicConstant = oray['AF_P'] / oray['PF_D']**oray['AF_Gm']

        fig2 = plt.figure()

        plt.semilogy( X1, PolytropicConstant, 'k-' )
        plt.axhline( self.EntropyThreshold, label = 'Shock radius cut-off' )

        plt.legend()
        plt.show()
        plt.close()

        del fig2

        input( '\nPress enter to proceed...' )

        return


    def MakeDataFile( self, OW_Option = False ):

        print( '\nCalling ShockRadius.MakeDataFile...\n' )

        if OW_Option:

            OW = True

        else:

            OW = OverwriteFile( self.DataFileName )

        if not OW: return

        FileArray        = self.FileArray
        nX               = self.nX
        DataDirectory    = self.DataDirectory
        PlotFileBaseName = self.PlotFileBaseName
        Field            = self.Field
        DataFileName     = self.DataFileName
        TimeFileName     = self.TimeFileName

        # Put all time-slices into one array to use for movie making

        Data = np.empty( (FileArray.shape[0],nX[0],nX[1]), np.float64 )

        Time = np.empty( (FileArray.shape[0]), np.float64 )

        print( '\nMaking data array for shock radius vs. time...' )
        print(   '----------------------------------------------' )

        for i in range( FileArray.shape[0] ):

            if (i+1) % 10 == 0:
                print( 'File {:}/{:}'.format( i+1, FileArray.shape[0] ) )

            Data[i], DataUnit, r, theta, Time[i], xL, xU \
              = GetData( DataDirectory, PlotFileBaseName, \
                         [ 'a', FileArray[i] ], Field, \
                         Verbose = False )

        np.savetxt( TimeFileName, Time )

        print( '\nMaking data file for shock radius vs. time...' )
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

                if i % 10 == 0:
                    print( 'File {:}/{:}'.format( i, FileArray.shape[0] ) )

                np.savetxt( FileOut, TimeSlice )
                FileOut.write( '# New slice\n' )

        return


    def ComputeShockRadius( self, Rs0 = 1.80e2 ):

        print( '\nCalling ShockRadius.ComputeShockRadius...\n' )

        OW = OverwriteFile( self.ShockRadiusVsTimeFileName )

        if not OW: return

        #self.MakeLineOutPlot()

        print( '\nComputing average shock radius...' )

        FileArray        = self.FileArray
        nX               = self.nX
        DataFileName     = self.DataFileName
        TimeFileName     = self.TimeFileName
        X1               = self.X1
        X2               = self.X2
        dX               = self.dX

        dX1 = dX[0]
        dX2 = dX[1]

        OW = OverwriteFile( DataFileName )

        self.MakeDataFile( OW_Option = OW )

        Data = np.loadtxt( DataFileName ).reshape( \
                 (FileArray.shape[0],nX[0],nX[1]) )
        Time = np.loadtxt( TimeFileName )

        DataUnit = 'erg/cm**3/(g/cm**3)**(4/3)'

        Volume = np.zeros( FileArray.shape[0], np.float64 )
        Min    = np.empty( FileArray.shape[0], np.float64 )
        Max    = np.empty( FileArray.shape[0], np.float64 )

        print( '\nComputing volume of shocked region...' )
        print(   '-------------------------------------' )

        for iT in range( Data.shape[0] ):

            if (iT+1) % 10 == 0:
                print( '  File {:}/{:}'.format( iT+1, Data.shape[0] ) )

            # Sum volumes of elements containing shocked fluid

            print( Data.min() )
            print( Data.max() );exit()
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

        print( '\nCalling ShockRadius.PlotShockRadiusVsTime...\n' )

        t, r, Min, Max = np.loadtxt( self.ShockRadiusVsTimeFileName )

        r0 = 1.0#[0]

        plt.plot( t, r  /r0, 'ko', label = 'Average' )
        plt.plot( t, Max/r0, 'bo', label = 'Max' )
        plt.plot( t, Min/r0, 'ro', label = 'Min' )
        #plt.xlim( t[0], t[-1] )
        plt.xlabel( 'Time [ms]' )
        plt.ylabel( r'$R_{s}/R_{s,0}$' )
        plt.legend()
        plt.show()
        plt.close()
        #plt.savefig( self.SaveFigAs, dpi = 300 )


if __name__ == "__main__":

    Root = '.'#'/Users/dunhamsj/Research/Data/AccretionShockParameterStudy/'
    EntropyThreshold = 1.0e15

    ID = 'NR2D_M1.4_Mdot0.3_Rs180_PA1.00e-04_nX640x064'
    Field = 'Entropy2'

    SR = ShockRadius( Root, ID, Field, EntropyThreshold )

    SR.ComputeShockRadius()
    SR.PlotShockRadiusVsTime()

    import os
    os.system( 'rm -rf __pycache__ ' )
