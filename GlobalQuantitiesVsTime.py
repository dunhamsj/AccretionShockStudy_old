#!/usr/bin/env python3

from scipy.optimize import curve_fit

import numpy as np
import matplotlib.pyplot as plt
import os
from sys import argv
plt.style.use( './Publication.sty' )

from UtilitiesModule import GetData as GD, GetFileArray, Overwrite

def FittingFunction( t, dlog10AKEdt, log10AKE0 ):

    return dlog10AKEdt * t + log10AKE0


class GlobalQuantities:

    def __init__( self, Root, ID ):

        self.Root = Root
        self.ID = ID

        self.DataDirectory = self.Root + self.ID + '/'
        self.PlotFileBaseName = self.ID + '.plt'

        return

    def IntegrateField( self, q, SqrtGm, dX1, dX2 ):

        nX1 = dX1.shape[0]
        nX2 = dX2.shape[0]

        Q = 0.0
        for iX1 in range( nX1 ):
            for iX2 in range( nX2 ):
                Q += 2.0 * np.pi * dX1[iX1] * dX2[iX2] \
                       * q[iX1,iX2] * SqrtGm[iX1,iX2]

        return Q

    def GetSnapshot( self, Field, iSS ):

        Data, DataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xU, nX, Time \
          = GD( self.DataDirectory, self.PlotFileBaseName, Field, \
                'spherical', True, argv = [ 'a', iSS ], \
                ReturnTime = True, ReturnMesh = True, Verbose = False )

        SqrtGm, DataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xU, nX, Time \
          = GD( self.DataDirectory, self.PlotFileBaseName, 'GF_SqrtGm', \
                'spherical', True, argv = [ 'a', iSS ], \
                ReturnTime = True, ReturnMesh = True, Verbose = False )

        dX1    *= 1.0e5
        SqrtGm *= 1.0e5**2

        return Data, SqrtGm, dX1, dX2, Time

    def GetData( self, Field ):

        FileName = self.ID + '_' + Field + '.dat'

        OW = Overwrite( FileName)#, ForceChoice = True, OW = False )

        if OW:

            FileArray \
              = GetFileArray( self.DataDirectory, self.PlotFileBaseName )#[::10]

            Q    = np.empty( FileArray.shape[0], np.float64 )
            Time = np.empty( FileArray.shape[0], np.float64 )

            for iSS in range( FileArray.shape[0] ):

                print( '{:d}/{:d}'.format( iSS+1, FileArray.shape[0] ) )

                Data, SqrtGm, dX1, dX2, Time[iSS] \
                  = self.GetSnapshot( Field, FileArray[iSS] )

                Q[iSS] = self.IntegrateField( Data, SqrtGm, dX1, dX2 )

            np.savetxt( FileName, np.vstack( ( Time, Q ) ) )

        d = np.loadtxt( FileName )

        Time = d[0]
        Q    = d[1]

        return Time, Q

    def FitCurve( self, Time, AKE ):

        beta, pcov = curve_fit( FittingFunction, Time, np.log10( AKE ) )

        return beta


if __name__ == '__main__':

    Root = '/lump/data/AccretionShockStudy/'
    IDs = np.array( [ 'NR2D_M2.0_Mdot0.3_Rs150', \
                      'GR2D_M2.0_Mdot0.3_Rs150' ], str )

    Fields = np.array( [ 'AngularKineticEnergyDensity' ], str )

    fig, axs = plt.subplots( Fields.shape[0], 1 )
    axs.set_title( IDs[0][5:] )
    axs.set_xlabel( 'Coordinate Time [ms]' )

    NR = GlobalQuantities( Root, IDs[0] )
    GR = GlobalQuantities( Root, IDs[1] )

    c = [ 'r', 'b' ]

    for i in range( Fields.shape[0] ):

        Time, dataNR = NR.GetData( Fields[i] )
        Time, dataGR = GR.GetData( Fields[i] )

        ind = np.where( Time > 10.0 )[0]

        Time   = np.copy( Time  [ind] )
        dataNR = np.copy( dataNR[ind] )
        dataGR = np.copy( dataGR[ind] )

        betaNR = NR.FitCurve( Time, dataNR )
        betaGR = GR.FitCurve( Time, dataGR )

        if i == 0:

            #axs.plot( Time, ( dataNR - dataNR[0] ) / dataNR[0], \
            #          c[0], label = 'NR' )
            #axs.plot( Time, ( dataGR - dataGR[0] ) / dataGR[0], \
            #          c[1], label = 'NR' )

            axs.plot( Time, dataNR, \
                      c[0], label = 'NR' )
            axs.plot( Time, dataGR, \
                      c[1], label = 'GR' )

            axs.plot( Time, 10**( betaNR[0] * Time + betaNR[1] ), c[0] + '--' )
            axs.plot( Time, 10**( betaGR[0] * Time + betaGR[1] ), c[1] + '--' )

            axs.text( 0.5, 0.3, 'm = {:.3e}'.format( betaNR[0] ), \
                      c = c[0], transform = axs.transAxes )
            axs.text( 0.5, 0.2, 'm = {:.3e}'.format( betaGR[0] ), \
                      c = c[1], transform = axs.transAxes )
            axs.text( 0.5, 0.1, r'(m$_{{GR}}$-m$_{{NR}}$)/m$_{{GR}}$ = {:.3e}' \
                      .format( ( betaGR[0] - betaNR[0] ) / betaGR[0] ), \
                      c = 'k', transform = axs.transAxes )

        else:

            #axs.plot( Time, ( dataNR - dataNR[0] ) / dataNR[0], \
            #          c[0] )
            #axs.plot( Time, ( dataGR - dataGR[0] ) / dataGR[0], \
            #          c[1] )

            axs.plot( Time, dataNR, \
                      c[0] )
            axs.plot( Time, dataGR, \
                      c[1] )

    axs.text( 0.2, 0.7, \
    r'$y\left(t\right)=2\pi\int\frac{W^{2}}{W+1}\,D\,v_{2}\,v^{2}\,\sqrt{\gamma}\,dr\,d\theta$', \
              transform = axs.transAxes )

    axs.set_ylabel( 'Angular Kinetic Energy erg/cm$^{3}$' )

    axs.legend()
    axs.set_yscale( 'log' )

    plt.savefig( 'fig.AKEComparison.png', dpi = 300 )
    #plt.show()

    os.system( 'rm -rf __pycache__' )
