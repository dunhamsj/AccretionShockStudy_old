#!/usr/bin/env python3

from scipy.integrate import simps
from scipy.optimize import curve_fit
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import os
from os.path import isfile
from sys import argv
plt.style.use( 'Publication.sty' )

from UtilitiesModule import ChoosePlotFile, Overwrite, GetFileArray

TwoPi = 2.0 * np.pi

class PowersInLegendreModes:

    def __init__( self, Root, ID, \
                  Field = 'Entropy', \
                  Rs = 1.80e2, fL = 0.90, fU = 0.95, R0 = -1.0, \
                  EntropyThreshold = 1.0e15, \
                  Verbose = True,
                  suffix = '', \
                  DataFileNameSuffix = '' ):

        self.Verbose = Verbose

        if self.Verbose:
            print( '\nCreating instance of PowersInLegendreModes class...\n' )

        self.Root  = Root
        self.ID    = ID
        self.Field = Field
        self.Rs    = Rs * 1.0e5
        self.fL    = fL
        self.fU    = fU
        self.R0    = R0 * 1.0e5
        self.EntropyThreshold = EntropyThreshold

        fMin = 40.0  / Rs
        fMax = 360.0 / Rs

        if self.fL < fMin:

            print( '\n  WARNING: fL < fMin. Setting fL to 40km / Rs\n' )

            self.fL = fMin

        if self.fU > fMax:

            print( '\n  WARNING: fU < fMax. Setting fU to 360km / Rs\n' )

            self.fU = fMax

        if self.R0 < 0.0:

            self.FileName \
              = '.{:}_{:}_{:.2f}-{:.2f}_PowersInLegendreModes{:}.dat'.format \
                  ( ID, Field, self.fL, self.fU, DataFileNameSuffix )

        else:

            self.FileName \
              = '.{:}_{:}_{:.2f}km_PowersInLegendreModes.dat'.format \
                  ( self.ID, self.Field, self.R0 / 1.0e5 )

        self.ShockRadiusVsTimeFileName \
          = '.{:}_ShockRadiusVsTime.dat'.format( self.ID )

        if self.Verbose:
            print( '  Variables:' )
            print( '  ----------' )
            print( '    Root:     {:s}'.format( self.Root ) )
            print( '    ID:       {:s}'.format( self.ID ) )
            print( '    Field:    {:s}'.format( self.Field ) )
            print( '    Rs:       {:.3e} km'.format( self.Rs / 1.0e5 ) )
            print( '    fL:       {:.2f}'.format( self.fL ) )
            print( '    fU:       {:.2f}'.format( self.fU ) )
            print( '    R0:       {:.3e} km'.format( self.R0 / 1.0e5 ) )
            print( '    Filename: {:s}\n'.format( self.FileName ) )

        self.suffix = suffix

        self.DataDirectory = self.Root + ID + '{:}/'.format( suffix )

        self.ComputedShockRadius = False
        self.ComputedPowers      = False

        return

    def FittingFunction( self, t, logF1, omega_r, omega_i, delta ):

        # Modified fitting function
        # (log of Eq. (9) in Blondin & Mezzacappa, (2006))

        return logF1 + 2.0 * omega_r * t \
                 + np.log( np.sin( omega_i * t + delta )**2 )

    def Jacobian( self, t, logF1, omega_r, omega_i, delta ):

        # Jacobian of modified fitting function

        J = np.empty( (t.shape[0],4), np.float64 )

        ImPhase = omega_i * t + delta

        J[:,0] = 1.0

        J[:,1] = 2.0 * t

        J[:,2] = 2.0 * np.cos( ImPhase ) / np.sin( ImPhase ) * t

        J[:,3] = 2.0 * np.cos( ImPhase ) / np.sin( ImPhase )

        return J

    def ComputeAngleAverage( self, Data, X2 ):
        return 1.0 / 2.0 * simps( Data * np.sin( X2 ), x = X2 )

    def GetShockRadiusVsTime( self ):

        Time, RsAve, RsMin, RsMax \
          = np.loadtxt( self.ShockRadiusVsTimeFileName )

        return Time, RsAve, RsMin, RsMax

    def ComputePowerInLegendreModes( self ):

        Time, RsAve, RsMin, RsMax, P0, P1, P2, P3, P4 \
          = np.loadtxt( self.FileName )

        return Time, RsAve, RsMin, RsMax, P0, P1, P2, P3, P4

    def FitPowerInLegendreModes \
      ( self, t, t0, t1, P1, \
        InitialGuess = np.array( [ 14.0, 1.0 / ( 2.0 * 200.0 ), TwoPi / 45.0, \
                                 0.0 ], np.float64 ) ):

        if self.Verbose:
            print( '\nCalling PowersInLegendreModes.FitPowerInLegendreModes...\n' )

        # --- Slice data for fitting ---

        ind = np.where( ( t >= t0 ) & ( t <= t1 ) )[0]

        tFit = t[ind] - t0

        # --- Fit model to data ---

        beta, pcov \
          = curve_fit( self.FittingFunction, tFit, np.log( P1[ind] ), \
                       p0 = InitialGuess, jac = self.Jacobian )
        F = np.exp( self.FittingFunction( tFit, beta[0], beta[1], \
                                          beta[2], beta[3] ) )

        self.beta = np.copy( beta )
        self.perr = np.sqrt( np.diag( pcov ) )

        # Propagate error from frequency into period
        self.perr[2] = TwoPi * self.perr[2] / self.beta[2]**2

        self.beta[0] = np.exp( self.beta[0] )
        self.beta[1] = self.beta[1] * 1.0e3
        self.perr[1] = self.perr[1] * 1.0e3
        self.beta[2] = TwoPi / self.beta[2]
        self.beta[3] = self.beta[3]

        return tFit + t0, F

    def PlotData( self, ax, m, rs, t0, t1, Time, P1, M, Rs, c, lw = 1 ):

        ind = np.where( ( Time >= t0 ) & ( Time <= t1 ) )[0]

        t = np.copy( Time )

        Time  = Time [ind]
        P1    = P1   [ind]

        if M > 0.0: ax.set_title( 'M{:}_Mdot0.3_Rs{:}'.format( M, Rs ) )

        if str( M ) == '2.8' and str( Rs ) == '120.0':
            if c == 'k-':
                ax.plot( Time, P1, c, lw = lw, label = 'WithPsi6' )
            elif c == 'r-':
                ax.plot( Time, P1, c, lw = lw, label = 'WithOutPsi6' )
            else:
                ax.plot( Time, P1, c, lw = lw, label = '|w/o - w|' )
        else:
            ax.plot( Time, P1, c, lw = lw )

        ax.set_yscale( 'log' )

        xlim = ( t0, t1 )
        ax.set_xlim( xlim )

        ax.set_xlabel( 'Time [ms]' )
        ax.set_ylabel( 'Power [cgs]' )

        ax.grid()

        return


if __name__ == "__main__":

    Root = '/lump/data/AccretionShockStudy/'

    Field = 'DivV2'
    t0    = 000.0
    t1    = 150.0
    fL    = 0.8
    fU    = 0.9
    R0    = -1.7e2
    suffix = ''

    M     = np.array( [ '1.4', '2.0', '2.8' ], str )
    Mdot  = np.array( [ '0.3' ], str )
    Rs    = np.array( [ '120', '150', '180' ], str )

    fig, axs = plt.subplots( 3, 3, figsize = (16,9) )
    fig.suptitle( Field + ' (GR)' )

    for mm in range( M.shape[0] ):
        for mdot in range( Mdot.shape[0] ):
            for rs in range( Rs.shape[0] ):

                m = M.shape[0] - mm - 1

                LogF0  = np.log( 1.0e14 )
                tauR   = 200.0
                T_SASI = 40.0
                delta  = 0.0

                tF0 = 1.0
                tF1 = 150.0

                if M[m] == '1.4':
                    if Rs[rs] == '120':
                        T_SASI = 25.0
                        tF0    = 15.0
                        tF1    = 120.0
                    elif Rs[rs] == '150':
                        T_SASI = 35.0
                        tF0    = 25.0
                        tF1    = 140.0
                    elif Rs[rs] == '180':
                        T_SASI = 55.0
                        tF0    = 35.0
                        tF1    = 140.0
                elif M[m] == '2.0':
                    if Rs[rs] == '120':
                        T_SASI = 20.0
                        tF0    = 1.0
                        tF1    = 150.0
                    elif Rs[rs] == '150':
                        T_SASI = 30.0
                        tF0    = 20.0
                        tF1    = 140.0
                    elif Rs[rs] == '180':
                        T_SASI = 50.0
                        tF0    = 1.0
                        tF1    = 150.0
                elif M[m] == '2.8':
                    if Rs[rs] == '120':
                        T_SASI = 20.0
                        tF0    = 15.0
                        tF1    = 150.0
                    elif Rs[rs] == '150':
                        T_SASI = 30.0
                        tF0    = 55.0
                        tF1    = 150.0
                    elif Rs[rs] == '180':
                        T_SASI = 40.0
                        tF0    = 5.0
                        tF1    = 150.0

                tF1 = 300.0
                omega_r = 1.0 / ( 2.0 * tauR )
                omega_i = TwoPi / T_SASI

                InitialGuess = np.array( [ LogF0, omega_r, omega_i, delta ], \
                                         np.float64 )

                ID_GR = 'GR2D_M{:}_Mdot{:}_Rs{:}'.format \
                        ( M[m], Mdot[mdot],  Rs[rs] )

                P_GR = PowersInLegendreModes( Root, ID_GR, Field, \
                                              Rs = np.float64( Rs[rs] ), \
                                              fL = fL, fU = fU, R0 = R0, \
                                              EntropyThreshold = 4.0e14, \
                                              Verbose = False )
                Time, RsAve, RsMin, RsMax, P0, P1, P2, P3, P4 \
                  = P_GR.ComputePowerInLegendreModes()
                P1_WithPsi6 = P1

                P_GR = PowersInLegendreModes( Root, ID_GR, Field, \
                                              Rs = np.float64( Rs[rs] ), \
                                              fL = fL, fU = fU, R0 = R0, \
                                              EntropyThreshold = 4.0e14, \
                                              Verbose = False, \
                                              DataFileNameSuffix = '_NoPsi6' )
                Time, RsAve, RsMin, RsMax, P0, P1, P2, P3, P4 \
                  = P_GR.ComputePowerInLegendreModes()

                P_GR.PlotData( axs[mm,rs], m, rs, t0, t1, Time, \
                               P1_WithPsi6, np.float64( M[m] ), \
                               np.float64( Rs[rs] ), 'k-', 2 )
                P_GR.PlotData( axs[mm,rs], m, rs, t0, t1, Time, \
                               P1         ,np.float64( M[m] ), \
                               np.float64( Rs[rs] ), 'r-' )
                P_GR.PlotData( axs[mm,rs], m, rs, t0, t1, Time, \
                               np.abs( P1 - P1_WithPsi6 ), np.float64( M[m] ), \
                               np.float64( Rs[rs] ), 'b-' )
                del ID_GR, P_GR, Time, RsAve, RsMin, RsMax, \
                    P0, P1, P2, P3, P4

    axs[0,0].legend()
    plt.subplots_adjust( wspace = 0.3, hspace = 0.4 )
    #plt.show()
    plt.savefig( \
    '/home/kkadoogan/fig.PowersInLegendreModes_MultiRun_{:}.png'.format \
    ( Field ), dpi = 300 )

    import os
    os.system( 'rm -rf __pycache__ ' )
