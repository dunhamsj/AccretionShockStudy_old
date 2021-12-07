#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import os

from MakeDataFile import MakeDataFile

class PlotFieldVsTheta:


    def __init__( self, Field, ID, DataDirectory, RunTime, nN, \
                  RadiiBelowShock, RadiiAboveShock ):

        self.Field           = Field
        self.ID              = ID
        self.DataDirectory   = DataDirectory
        self.MovieRunTime    = RunTime
        self.nN              = nN
        self.RadiiBelowShock = RadiiBelowShock
        self.RadiiAboveShock = RadiiAboveShock

        self.SaveFileAs_Nodal = 'mov.{:}_{:}_Nodal.mp4'.format( ID, Field )

        return


    def GetDomainInfo( self ):

        ID = self.ID
        Field = self.Field
        DataDirectory = self.DataDirectory + ID + '/'
        nN = self.nN

        PlotFileBaseName = ID + '.plt'

        self.DataFileName = '.' + ID + '_' + Field + '.dat'

        PlotFileBaseName = ID + '.plt'

        xL, xU, nX, FileArray \
          = MakeDataFile( Field, DataDirectory, \
                          self.DataFileName, PlotFileBaseName )

        self.xL = xL
        self.xU = xU
        self.nX = nX
        self.FileArray = FileArray

        nSS = FileArray.shape[0]
        self.nSS = nSS

        dX = ( xU - xL ) / np.float64( nX )
        self.dX = dX

        nX[0] += 2
        nX[1] += 2

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

        # Create theta mesh

        thetaC = np.linspace( xL[1] - 0.5 * dX[1], xU[1] + 0.5 * dX[1], nX[1] )
        theta  = np.empty( (nDOFX[1]), np.float64 )

        for i in range( nX[1] ):
            ind = np.linspace( nN*i, nN*i+nN-1, nN, dtype = np.int64 )
            theta[ind] = xq * dX[1] + thetaC[i]

        self.theta = theta

        # Create radial mesh

        rC = np.linspace( xL[0] - 0.5 * dX[0], xU[0] + 0.5 * dX[0], nX[0] )
        r  = np.empty( (nDOFX[0]), np.float64 )

        for i in range( nX[0] ):
            ind = np.linspace( nN*i, nN*i+nN-1, nN, dtype = np.int64 )
            r[ind] = xq * dX[0] + rC[i]

        self.r = r

        return


    def GetData_Nodal( self ):

        DataFileName = self.DataFileName
        RadiiBelowShock = self.RadiiBelowShock
        RadiiAboveShock = self.RadiiAboveShock
        dX = self.dX
        r = self.r
        nX = self.nX
        nN = self.nN
        nSS = self.nSS
        FileArray = self.FileArray
        DataDirectory = self.DataDirectory + self.ID + '/'

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

        ind0 = np.empty( (RadiiBelowShock.shape[0]), np.int64 )
        ind1 = np.empty( (RadiiAboveShock.shape[0]), np.int64 )

        for i in range( ind0.shape[0] ):
            ind0[i] \
              = np.where(  ( r > RadiiBelowShock[i] - 0.5 * dX[0] ) \
                         & ( r < RadiiBelowShock[i] + 0.5 * dX[0] ) )[0][0]

        for i in range( ind1.shape[0] ):
            ind1[i] \
              = np.where(   ( r > RadiiAboveShock[i] - 0.5 * dX[0] ) \
                          & ( r < RadiiAboveShock[i] + 0.5 * dX[0] ) )[0][0]

        self.ind0 = ind0
        self.ind1 = ind1

        nX0 = np.int64( ind0 / nN )
        nX1 = np.int64( ind1 / nN )

        indT0 = ind0 % nN
        indT1 = ind1 % nN

        d = np.empty( (nX[0],nX[1],nN,nN), np.float64 )

        V20 = np.empty( (ind0.shape[0],nSS,nX[1]*nN), np.float64 )
        V21 = np.empty( (ind1.shape[0],nSS,nX[1]*nN), np.float64 )

        for i in range( nSS ):

            if (i+1) % 10 == 0:
                print( 'File {:d}/{:d}'.format( i+1, nSS ) )

            FileNumber = FileArray[i][-8:]

            Data = np.loadtxt( DataDirectory + 'NodalData_TCI0_{:}_V2.dat'.format \
                               ( FileNumber ) )

            j = 0
            for iX2 in range( nX[1] ):
                for iX1 in range( nX[0] ):

                    if   iX1 == 0       and iX2 == 0:       d[iX1,iX2] = 0.0
                    elif iX1 == 0       and iX2 == nX[1]-1: d[iX1,iX2] = 0.0
                    elif iX1 == nX[0]-1 and iX2 == 0:       d[iX1,iX2] = 0.0
                    elif iX1 == nX[0]-1 and iX2 == nX[1]-1: d[iX1,iX2] = 0.0
                    else:

                        d[iX1,iX2] = Data[j].reshape( (nN,nN) )
                        j += 1

            for j in range( ind0.shape[0] ):

                V20[j,i] = d[nX0[j],:,:,indT0[j]].flatten()

            for j in range( ind1.shape[0] ):

                V21[j,i] = d[nX1[j],:,:,indT1[j]].flatten()

        self.V20_Nodal = V20
        self.V21_Nodal = V21

        return


    def MakeMovie1D( self ):

        xL = self.xL
        xU = self.xU
        dX = self.dX
        V20 = self.V20_Nodal
        V21 = self.V21_Nodal
        Time = self.Time
        Field = self.Field
        FieldUnit = self.FieldUnit
        ind0 = self.ind0
        ind1 = self.ind1
        r = self.r
        theta = self.theta
        nSS = self.nSS
        MovieRunTime = self.MovieRunTime
        SaveFileAs = self.SaveFileAs_Nodal

        fig, axs = plt.subplots( 2, 1 )

        xlim  = [ xL[1]-2.0*dX[1], xU[1]+2.0*dX[1] ]
        Width = xlim[1] - xlim[0]

        ylim0 = [ V20.min(), V20.max() ]
        ylim1 = [ V21.min(), V21.max() ]
        Height = ylim0[1] - ylim0[0]

        axs[0].set_xlim( xlim )
        axs[1].set_xlim( xlim )

        axs[0].set_ylim( ylim0 )
        axs[1].set_ylim( ylim1 )

        axs[0].set_xlabel( r'$\theta$' )
        axs[1].set_xlabel( r'$\theta$' )

        axs[0].set_ylabel( '{:} {:}'.format( Field, FieldUnit ) )
        axs[1].set_ylabel( '{:} {:}'.format( Field, FieldUnit ) )

        time_text = axs[0].text( xlim [0] + 0.6 * Width, \
                                 ylim0[0] + 0.9 * Height, '' )

        lines0 = np.empty( ind0.shape[0], object )
        lines1 = np.empty( ind1.shape[0], object )

        for i in range( lines0.shape[0] ):
            lines0[i], = axs[0].plot( [], [], 'k.', markersize = 1.0, \
                              label = 'r = {:.1f} km'.format( r[ind0[i]] ) )

        for i in range( lines1.shape[0] ):
            lines1[i], = axs[1].plot( [], [], 'k.', markersize = 1.0, \
                              label = 'r = {:.1f} km'.format( r[ind1[i]] ) )

        axs[0].legend( loc = 2 )
        axs[1].legend( loc = 2 )

        def InitializeFrame():

            ret = []

            for i in range( lines0.shape[0] ):
                lines0[i].set_data([],[])
                ret.append( lines0[i] )

            for i in range( lines1.shape[0] ):
                lines1[i].set_data([],[])
                ret.append( lines1[i] )

            time_text.set_text('')
            ret.append( time_text )
            ret = tuple( ret )

            return ret

        def UpdateFrame(t):

            ret = []

            for i in range( lines0.shape[0] ):
                y = V20[i,t]
                lines0[i].set_data( theta, y )
                ret.append( lines0[i] )

            for i in range( lines1.shape[0] ):
                y = V21[i,t]
                lines1[i].set_data( theta, y )
                ret.append( lines1[i] )

            time_text.set_text( 'time = {:d} ms'.format( np.int64( Time[t] ) ) )
            ret.append( time_text )

            ret = tuple( ret )
            return ret

        fps = np.float64( nSS ) / MovieRunTime

        plt.subplots_adjust( hspace = 0.4 )

        anim = animation.FuncAnimation( fig, UpdateFrame, \
                                        init_func = InitializeFrame, \
                                        frames = nSS, \
                                        blit = True )

        anim.save( SaveFileAs, fps = fps, dpi = 300 )

        os.system( 'rm -rf __pycache__ ' )

        return

if __name__ == "__main__":

    Field = 'PF_V2'

    ID = 'GR2D_M1.4_Mdot0.3_Rs180_PA0.000_nX320x064_TCI0'

    DataDirectory \
      = '/Users/dunhamsj/Research/Data/AccretionShockParameterStudy/'

    RunTime = 10.0

    nN = 3

    RadiiBelowShock = np.array( [ 100.0 ], np.float64 )
    RadiiAboveShock = np.array( [ 300.0 ], np.float64 )

    P = PlotFieldVsTheta( Field, ID, DataDirectory, RunTime, nN, \
                          RadiiBelowShock, RadiiAboveShock )

    P.GetDomainInfo()
    P.GetData_Nodal()
    P.MakeMovie1D()

