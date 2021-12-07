#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import os

from MakeDataFile import MakeDataFile

class NodalData:

    DataDirectory \
      = '/Users/dunhamsj/Research/Data/AccretionShockParameterStudy/'

    def __init__( self, Field, UseLogScale, ID, MovieRunTime, \
                  nN, ShockRadius, nBorder, SSi, SSf ):

        self.Field         = Field
        self.UseLogScale   = UseLogScale
        self.ID            = ID
        self.MovieRunTime  = MovieRunTime
        self.nN            = nN
        self.ShockRadius   = ShockRadius
        self.nBorder       = nBorder
        self.SSi           = SSi
        self.SSf           = SSf

        self.DataDirectory += ID + '_DS1E/'
        self.SaveFileAs     = 'mov.{:}_{:}_Nodal.mp4'.format( ID, Field )

        return


    def GetDomainInfo( self ):

        ID = self.ID
        Field = self.Field
        nN = self.nN
        SSi = self.SSi
        SSf = self.SSf

        PlotFileBaseName = ID + '.plt'

        self.DataFileName = '.' + ID + '_' + Field + '.dat'

        xL, xU, nX, FileArray \
          = MakeDataFile( Field, self.DataDirectory, \
                          self.DataFileName, PlotFileBaseName, \
                          SSi = self.SSi, SSf = self.SSf )

        self.xL = xL
        self.xU = xU
        self.nX = nX
        self.FileArray = FileArray

        nSS = FileArray.shape[0]-1
        self.nSS = nSS

        dX = ( xU - xL ) / np.float64( nX )
        self.dX = dX

        nX[0] += 2

        self.nX = nX

        nDOFX = nX * nN
        self.nDOFX = nDOFX

        if   nN == 2:

            xq = np.array( [ -np.sqrt( 1.0 / 12.0 ), \
                             +np.sqrt( 1.0 / 12.0 ) ], \
                           np.float64 )

        elif nN == 3:

            xq = np.array( [ -np.sqrt( 3.0 / 20.0 ), 0.0, \
                             +np.sqrt( 3.0 / 20.0 ) ], \
                           np.float64 )

        # Create radial mesh

        rC = np.linspace( xL[0] - 0.5 * dX[0], xU[0] + 0.5 * dX[0], nX[0] )
        r  = np.empty( (nDOFX[0]), np.float64 )

        for i in range( nX[0] ):
            ind = np.linspace( nN*i, nN*i+nN-1, nN, dtype = np.int64 )
            r[ind] = xq * dX[0] + rC[i]

        self.r = r
        self.rC = rC

        return


    def GetData_Nodal( self, SSi = -1, SSf = -1 ):

        # Get header data

        if SSi < 0: SSi = 0
        if SSf < 0: SSf = self.nSS

        self.nSS = SSf - SSi

        Field = self.Field

        DataFileName = self.DataFileName
        nBorder = self.nBorder
        dX = self.dX
        r  = self.r
        rC = self.rC
        nX = self.nX
        nN = self.nN
        nSS = self.nSS
        FileArray = self.FileArray
        ShockRadius = self.ShockRadius

        f = open( DataFileName )
        header = f.readline()[16:-2]
        FieldUnit = f.readline()[9:-1]
        self.FieldUnit = FieldUnit
        DataShape = tuple( [ np.int64( dim ) for dim in header.split( ',' ) ] )

        Time = list( [ np.float64( t ) \
                       for t in f.readline()[12:-2].split(' ') ] )
        Time = np.array( Time )

        self.Time = Time

        Data = np.loadtxt( DataFileName, dtype = np.float64 ).reshape \
                ( DataShape )

        # Find radial DOF corresponding to element

        rB = ShockRadius - nBorder * dX[0]
        rA = ShockRadius + nBorder * dX[0]

        ind      = np.where( ( r > rB ) & ( r < rA ) )[0]
        self.ind = ind

        Data = np.empty( (nSS,ind.shape[0]), np.float64 )

        for i in range( SSi, SSf ):

            j = i - SSi

            if (j+1) % 10 == 0:
                print( 'File {:d}/{:d}'.format( j+1, nSS ) )

            FileNumber = FileArray[i][-8:]

            d = np.loadtxt( self.DataDirectory \
                            + ID + '.Nodal.{:}.{:}.dat'.format \
                            ( FileNumber, Field ) ).flatten()

            Data[j] = d[ind]

        self.Data = Data

        return


    def MakeMovie1D( self ):

        xL = self.xL
        xU = self.xU
        dX = self.dX
        Data = self.Data
        Time = self.Time
        Field = self.Field
        FieldUnit = self.FieldUnit
        ind = self.ind
        r = self.r
        nSS = self.nSS-1

        MovieRunTime = self.MovieRunTime
        SaveFileAs = self.SaveFileAs

        fig, ax = plt.subplots( 1, 1, figsize = (8,6) )

        fig.suptitle( Field, y = 0.95 )

        rN = self.r[self.ind]

        xlim  = [ rN[0]-2.0*dX[0], rN[-1]+2.0*dX[0] ]
        Width = xlim[1] - xlim[0]

        D    = np.copy( Data )
        Data = np.empty( (D.shape[0]-1,D.shape[1]), np.float64 )

        for iSS in range( nSS ):

            for iNX in range( D.shape[1] ):

                Data[iSS,iNX] \
                  = max( 1.0e-17, np.abs( ( D[iSS+1,iNX] - D[iSS,iNX] ) \
                          / ( 0.5 * ( D[iSS+1,iNX] + D[iSS,iNX] ) ) ) )

        ylim = [ Data.min(), Data.max() ]
        Height = ylim[1] - ylim[0]

        ax.set_xlim( xlim )
        ax.set_ylim( ylim )

        ax.set_xlabel( 'r km' )

        ax.grid()

        ax.set_ylabel( r'$\left|\frac{y^{n+1}-y^{n}}{\frac{1}{2}\left(y^{n+1}+y^{n}\right)}\right|$', fontsize = 15 )

        if self.UseLogScale: ax.set_yscale( 'log' )

        time_text = ax.text( 0.8, 1.01, '', transform=ax.transAxes )

        lines, = ax.plot( [], [], 'k.', markersize = 1.0 )

        def InitializeFrame():

            ret = []

            lines.set_data([],[])
            ret.append( lines )

            time_text.set_text('')
            ret.append( time_text )

            return ( ret )

        def UpdateFrame(t):

            ret = []

            lines.set_data( rN, Data[t] )
            ret.append( lines )

            time_text.set_text( 'time = {:d} ms'.format( np.int64( Time[t] ) ) )
            ret.append( time_text )

            return ( ret )

        fps = np.float64( nSS ) / MovieRunTime

        plt.subplots_adjust( hspace = 0.4 )

        anim = animation.FuncAnimation( fig, UpdateFrame, \
                                        init_func = InitializeFrame, \
                                        frames = nSS, \
                                        blit = True )

        anim.save( SaveFileAs, fps = fps, dpi = 300 )

        os.system( 'rm -rf __pycache__ ' )

        return


    def MakeMovie1D_Conserved( self, D, S, E, PlotDifference = False ):

        xL   = self.xL
        xU   = self.xU
        dX   = self.dX
        Time = self.Time
        ind  = self.ind
        r    = self.r
        nSS  = self.nSS-1

        MovieRunTime = self.MovieRunTime

        fig, axs = plt.subplots( 3, 1, figsize = (8,6) )

#        fig.suptitle( Field, y = 0.95 )

        rN = self.r[self.ind]

        xlim  = [ rN[0]-2.0*dX[0], rN[-1]+2.0*dX[0] ]
        Width = xlim[1] - xlim[0]

        axs[0].set_ylabel( r'$D\,\mathrm{g/cm}^{3}$' )
        axs[1].set_ylabel( r'$S_{1}\,\mathrm{g/s/cm}^{2}$' )
        axs[2].set_ylabel( r'$E\,\mathrm{erg/cm}^{3}$' )

        if PlotDifference:

            yD = np.copy( D )
            yS = np.copy( S )
            yE = np.copy( E )

            D = np.empty( (D.shape[0]-1,D.shape[1]), np.float64 )
            S = np.empty( (S.shape[0]-1,S.shape[1]), np.float64 )
            E = np.empty( (E.shape[0]-1,E.shape[1]), np.float64 )

            nDOFX = D.shape[1]

            for iSS in range( nSS ):

                for iNX in range( nDOFX ):

                    D[iSS,iNX] \
                      = max( 1.0e-17, np.abs( ( yD[iSS+1,iNX] - yD[iSS,iNX] ) \
                              / ( 0.5 * ( yD[iSS+1,iNX] + yD[iSS,iNX] ) ) ) )

                    S[iSS,iNX] \
                      = max( 1.0e-17, np.abs( ( yS[iSS+1,iNX] - yS[iSS,iNX] ) \
                              / ( 0.5 * ( yS[iSS+1,iNX] + yS[iSS,iNX] ) ) ) )

                    E[iSS,iNX] \
                      = max( 1.0e-17, np.abs( ( yE[iSS+1,iNX] - yE[iSS,iNX] ) \
                              / ( 0.5 * ( yE[iSS+1,iNX] + yE[iSS,iNX] ) ) ) )

            axs[0].set_ylabel( r'$\left|\frac{D^{n+1}-D^{n}}{\frac{1}{2}\left(D^{n+1}+D^{n}\right)}\right|$', fontsize = 15 )

            axs[1].set_ylabel( r'$\left|\frac{S^{n+1}-S^{n}}{\frac{1}{2}\left(S^{n+1}+S^{n}\right)}\right|$', fontsize = 15 )

            axs[2].set_ylabel( r'$\left|\frac{E^{n+1}-E^{n}}{\frac{1}{2}\left(E^{n+1}+E^{n}\right)}\right|$', fontsize = 15 )

        ylim = [ [ D.min(), D.max() ], \
                 [ S.min(), S.max() ], \
                 [ E.min(), E.max() ] ]

        Height = []

        for i in range( 3 ):

            Height.append( ylim[i][1] - ylim[i][0] )
            axs[i].set_xlim( xlim )
            axs[i].set_ylim( ylim[i] )
            if self.UseLogScale[i]: axs[i].set_yscale( 'log' )
            axs[i].grid( axis = 'x' )

        axs[2].set_xlabel( 'r km' )

        time_text = axs[0].text( 0.8, 1.01, '', transform = axs[0].transAxes )
        lineD, = axs[0].plot( [], [], 'k.', markersize = 1.0 )
        lineS, = axs[1].plot( [], [], 'k.', markersize = 1.0 )
        lineE, = axs[2].plot( [], [], 'k.', markersize = 1.0 )

        def InitializeFrame():

            ret = []

            lineD.set_data([],[])
            lineS.set_data([],[])
            lineE.set_data([],[])
            ret.append( lineD )
            ret.append( lineS )
            ret.append( lineE )

            time_text.set_text('')
            ret.append( time_text )

            return ( ret )

        def UpdateFrame(t):

            ret = []

            lineD.set_data( rN, D[t] )
            lineS.set_data( rN, S[t] )
            lineE.set_data( rN, E[t] )
            ret.append( lineD )
            ret.append( lineS )
            ret.append( lineE )

            time_text.set_text( 'time = {:d} ms'.format( np.int64( Time[t] ) ) )
            ret.append( time_text )

            return ( ret )

        fps = np.float64( nSS ) / MovieRunTime

        plt.subplots_adjust( hspace = 0.0 )

        anim = animation.FuncAnimation( fig, UpdateFrame, \
                                        init_func = InitializeFrame, \
                                        frames = nSS, \
                                        blit = True )

        SaveFileAs = 'mov.MaxRelDiff.Conserved.mp4'
        anim.save( SaveFileAs, fps = fps, dpi = 300 )

        os.system( 'rm -rf __pycache__ ' )

        return

if __name__ == "__main__":

    UseLogScale  = [ False, False, False ]
    ID           = 'GR1D_M1.4_Mdot0.3_Rs180_PA0.000_nX320'
    MovieRunTime = 10.0
    nN           = 3
    ShockRadius  = 1.80e2
    nBorder      = 10

    SSi = -1
    SSf = -1

    P_D = NodalData( 'CF_D', UseLogScale, ID, \
                     MovieRunTime, nN, ShockRadius, nBorder, SSi, SSf )
    P_D.GetDomainInfo()
    P_D.GetData_Nodal( SSi, SSf )

    P_S = NodalData( 'CF_S1', UseLogScale, ID, \
                     MovieRunTime, nN, ShockRadius, nBorder, SSi, SSf )
    P_S.GetDomainInfo()
    P_S.GetData_Nodal( SSi, SSf )

    P_E = NodalData( 'CF_E', UseLogScale, ID, \
                     MovieRunTime, nN, ShockRadius, nBorder, SSi, SSf )
    P_E.GetDomainInfo()
    P_E.GetData_Nodal( SSi, SSf )

    P_D.MakeMovie1D_Conserved( P_D.Data, P_S.Data, P_E.Data )

#    P_D.MakeMovie1D()

