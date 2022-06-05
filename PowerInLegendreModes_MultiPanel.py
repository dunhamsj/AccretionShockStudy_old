#!/usr/bin/env python3

import yt
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

yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen
TwoPi = 2.0 * np.pi

def FittingFunction2( t, logF1, omega_r, omega_i, delta ):

    # Modified fitting function
    # (log of Eq. (9) in Blondin & Mezzacappa, (2006))

    return logF1 + 2.0 * omega_r * t \
             + np.log( np.sin( omega_i * t + delta )**2 )

def Jacobian2( t, logF1, omega_r, omega_i, delta ):

    # Jacobian of modified fitting function

    J = np.empty( (t.shape[0],4), np.float64 )

    ImPhase = omega_i * t + delta

    J[:,0] = 1.0

    J[:,1] = 2.0 * t

    J[:,2] = 2.0 * np.cos( ImPhase ) / np.sin( ImPhase ) * t

    J[:,3] = 2.0 * np.cos( ImPhase ) / np.sin( ImPhase )

    return J

class PowersInLegendreModes:

    def __init__( self, Root, ID, \
                  Field = 'Entropy', \
                  Rs = 1.80e2, fL = 0.90, fU = 0.95, R0 = -1.0, \
                  EntropyThreshold = 1.0e15, \
                  Verbose = True,
                  suffix = '' ):

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
              = '.{:}_{:}_{:.2f}-{:.2f}_PowersInLegendreModes.dat'.format \
                  ( ID, Field, self.fL, self.fU )

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

    def FittingFunction( self, t, beta ):

        # Modified fitting function
        # (log of Eq. (9) in Blondin & Mezzacappa, (2006))

        logF1   = beta[0]
        omega_r = beta[1]
        omega_i = beta[2]
        delta   = beta[3]

        return logF1 + 2.0 * omega_r * t \
                 + np.log( np.sin( omega_i * t + delta )**2 )

    def Jacobian( self, t, beta ):

        # Jacobian of modified fitting function

        J = np.empty( (t.shape[0],4), np.float64 )

        logF1   = beta[0]
        omega_r = beta[1]
        omega_i = beta[2]
        delta   = beta[3]

        ImPhase = omega_i * t + delta

        J[:,0] = 1.0

        J[:,1] = 2.0 * t

        J[:,2] = 2.0 * np.cos( ImPhase ) / np.sin( ImPhase ) * t

        J[:,3] = 2.0 * np.cos( ImPhase ) / np.sin( ImPhase )

        return J

    def CurveFit( self, t, logP, InitialGuess ):

        if self.Verbose:
            print( '\nCalling PowersInLegendreModes.CurveFit...\n' )

        nIter = 1000

        beta_k = np.copy( InitialGuess )

        eps  = 1.0e-10
        db   = 1.0
        iter = 0

        alpha = 0.5

        while( db > eps and iter < nIter ):

            iter += 1

            F = self.FittingFunction( t, beta_k )
            J = self.Jacobian       ( t, beta_k )

            r = logP - F # Residual

            dbeta = alpha * np.dot( inv( np.dot( J.T, J ) ), \
                            np.dot( J.T, r ) )

            beta_k += dbeta

            db = np.max( np.abs( dbeta / beta_k ) )

#            if iter > nIter - 5:
#
#                print( 'iter = {:d}, db = {:.3e}'.format( iter, db ) )

        print( 'iter = {:d}, db = {:.3e}'.format( iter, db ) )

        return beta_k

    def GetPowersFromDiagnostics( self ):

        if self.Verbose:
            print( '\nCalling PowersInLegendreModes.GetPowersFromDiagnostics...\n' )

        ID = self.ID

        # --- Read in diagnostics ---

        Data = np.loadtxt( self.DataDirectory + ID + '.Diagnostics.dat', \
                           skiprows = 1 )

        t     = np.copy( Data[:,0] )
        P0    = np.copy( Data[:,1] )
        P1    = np.copy( Data[:,2] )
        P2    = np.copy( Data[:,3] )
        RsAve = np.copy( Data[:,4] )
        RsMin = np.copy( Data[:,5] )
        RsMax = np.copy( Data[:,6] )

        return t, P0, P1, P2, RsAve, RsMin, RsMax

    def ComputeAngleAverage( self, Data, X2 ):
        return 1.0 / 2.0 * simps( Data * np.sin( X2 ), x = X2 )

    def GetShockRadiusVsTime( self ):

        if not self.ComputedShockRadius:

            if self.Verbose:
                print( '\nCalling PowersInLegendreModes.GetShockRadiusVsTime...\n' )

            OW = Overwrite( self.ShockRadiusVsTimeFileName, \
                            ForceChoice = True, OW = False )

            if OW:

                from ShockRadius import ShockRadius
                SR = ShockRadius( self.Root, self.ID, \
                                  EntropyThreshold = self.EntropyThreshold, \
                                  suffix = self.suffix )
                SR.ComputeShockRadius()

            self.ComputedShockRadius = True

        Time, RsAve, RsMin, RsMax \
          = np.loadtxt( self.ShockRadiusVsTimeFileName )

        return Time, RsAve, RsMin, RsMax

    def ComputePowerInLegendreModes( self ):

        Field = self.Field

        ID = self.ID

        PlotFileBaseName = ID + '.plt'

        FileArray = GetFileArray( self.DataDirectory, PlotFileBaseName )

        nFiles = FileArray.shape[0]

        # Get info about computational domain

        ds = yt.load( '{:}'.format( self.DataDirectory + FileArray[0] ) )

        MaxLevel = ds.index.max_level
        nX       = ds.domain_dimensions
        xL       = ds.domain_left_edge.to_ndarray()
        xH       = ds.domain_right_edge.to_ndarray()

        CentimetersPerKilometer = 1.0e5

        xL[0] *= CentimetersPerKilometer
        xH[0] *= CentimetersPerKilometer

        dX = ( xH - xL ) / np.float64( nX )

        X1 = np.linspace( xL[0] + dX[0] / 2.0, xH[0] - dX[0] / 2.0, nX[0] )
        X2 = np.linspace( xL[1] + dX[1] / 2.0, xH[1] - dX[1] / 2.0, nX[1] )

        x = np.cos( X2 )

        # Legendre polynomials

        P0 = np.sqrt( 1.0 / 2.0 ) \
               * np.ones( nX[1] )
        P1 = np.sqrt( 3.0 / 2.0 ) \
               * x
        P2 = np.sqrt( 5.0 / 2.0 ) \
               * ( 3.0 * x**2 - 1.0 ) / 2.0
        P3 = np.sqrt( 7.0 / 2.0 ) \
               * 1.0 / 2.0 * ( 5.0 * x**3 - 3.0 * x )
        P4 = np.sqrt( 9.0 / 2.0 ) \
               * 1.0 / 8.0 * ( 35.0 * x**4 - 30.0 * x**2 + 3.0 )

        # Radially-dependent quantity (G)

        G0 = np.empty( (nFiles,nX[0]), np.float64 )
        G1 = np.empty( (nFiles,nX[0]), np.float64 )
        G2 = np.empty( (nFiles,nX[0]), np.float64 )
        G3 = np.empty( (nFiles,nX[0]), np.float64 )
        G4 = np.empty( (nFiles,nX[0]), np.float64 )

        # Powers in Legendre modes

        Power = np.empty( (nFiles,5), np.float64 )

        if self.ComputedPowers:

            Time, RsAve, RsMin, RsMax, P0, P1, P2, P3, P4 \
              = np.loadtxt( self.FileName )

            return Time, RsAve, RsMin, RsMax, P0, P1, P2, P3, P4

        if self.Verbose:
            print( \
              '\nCalling PowersInLegendreModes.ComputePowerInLegendreModes...' )
            print( \
              '------------------------------------------------------------\n' )

        OW = Overwrite( self.FileName)#, ForceChoice = True, OW = True )

        self.ComputedPowers = True

        if OW:

            Time = np.empty( nFiles )

            for i in range( nFiles ):

                if (i+1) % 10 == 0:
                    print( 'File {:}/{:}'.format( i+1, nFiles ) )

                ds = yt.load( '{:}'.format( self.DataDirectory \
                                              + FileArray[i] ) )
                Time[i] = ds.current_time

                CoveringGrid \
                  = ds.covering_grid \
                      ( level           = ds.index.max_level, \
                        left_edge       = ds.domain_left_edge, \
                        dims            = ds.domain_dimensions * 2**MaxLevel, \
                        num_ghost_zones = ds.domain_dimensions[0] )

                ds.force_periodicity()

                PF_V1 = CoveringGrid['PF_V1'].to_ndarray()[:,:,0] \
                          * CentimetersPerKilometer
                PF_V2 = CoveringGrid['PF_V2'].to_ndarray()[:,:,0]
                PF_D  = CoveringGrid['PF_D' ].to_ndarray()[:,:,0]
                AF_P  = CoveringGrid['AF_P' ].to_ndarray()[:,:,0]
                AF_Gm = CoveringGrid['AF_Gm'].to_ndarray()[:,:,0]

                fL = self.fL
                fU = self.fU
                Rs = self.Rs

                if self.R0 < 0.0:

                    indX1 \
                      = np.where( ( X1 > fL * Rs ) & ( X1 < fU * Rs ) )[0]

                else:

                    indX1 \
                      = np.where(   ( X1 > self.R0 - dX[0] ) \
                                  & ( X1 < self.R0 + dX[0]) )[0]

                    indX1 = np.array( [ indX1[0] ], np.int64 )

                indX2 = np.linspace( 0, nX[1]-1, nX[1], dtype = np.int64 )

                Data = np.zeros( (nX[0],nX[1]), np.float64 )

                if  ( Field == 'OverPressure' ):

                    AngleAveragedPressure \
                      = np.zeros( (nX[0],nX[1]), np.float64 )

                    for j in indX1:
                        AngleAveragedPressure[j]\
                          = ComputeAngleAverage( AF_P[j,:], X2 )

                    # Loop over radius

                    for j in indX1:
                        Data[j,:] = ( AF_P[j,:] - AngleAveragedPressure[j] ) \
                                      / AngleAveragedPressure[j]

                elif( Field == 'RatioOfKineticEnergies' ):

                    for j in indX1:
                        Data[j,:] = ( X1[j] * PF_V2[j,:] / PF_V1[j,:] )**2

                elif( Field == 'Entropy' ):

                    Data = np.log10( AF_P / PF_D**(AF_Gm) )

                elif( Field == 'V2' ):

                    Data = PF_V2

                elif( Field == 'pr4' ):

                    for j in indX1:
                        for k in indX2:

                            Data[j,k] = AF_P[j,k] * X1[j]**4

                elif( Field == 'DivV2' ):

                    # --- Sheck et al., (2008), A&A, 477, 931 ---

                    indX2 = np.linspace( 1, nX[1]-2, nX[1]-2, dtype = np.int64 )

                    for j in indX1:
                        for k in indX2:

                            Data[j,k] \
                              = 1.0 / ( 2.0 * dX[1] * np.sin( X2[k] ) ) \
                                  * (   np.sin( X2[k+1] ) * PF_V2[j,k+1] \
                                      - np.sin( X2[k-1] ) * PF_V2[j,k-1] )

                elif( Field == 'Vorticity' ):

                    indX1 = np.linspace( 0, nX[0]-2, nX[0]-1, dtype = np.int64 )
                    indX2 = np.linspace( 0, nX[1]-3, nX[1]-2, dtype = np.int64 )

                    for j in indX1:
                        for k in indX2:

                            Data[j,k] \
                              = 1.0 / ( 1.0 / 2.0 * ( X1[j+1] + X1[j] ) ) \
                                  * ( ( X1[j+1]**2 * PF_V2[j+1,k+2] \
                                          - X1[j]**2 * PF_V2[j,k+1] ) / dX1 \
                                       - ( PF_V1[j,k+2] - PF_V1[j,k+1] ) / dX2 )

                else:

                    print( 'Invalid choice of field: {:s}'.format( Field ) )
                    print( 'Valid choices' )
                    print( '-------------' )
                    print( '  OverPressure' )
                    print( '  RatioOfKineticEnergies' )
                    print( '  Entropy' )
                    print( '  V2' )
                    print( '  pr4' )
                    print( '  DivV2' )
                    print( '  Vorticity' )
                    exit( 'Exiting...' )

#                # Subtract off angle-average
#                for j in indX1:
#                    Data[j,indX2] \
#                      -= self.ComputeAngleAverage( Data[j,indX2], X2[indX2] )

                # --- For each radius, j, integrate over full theta range ---

                for j in indX1:

                    G0[i,j] \
                      = simps( Data[j,indX2] * P0[indX2] \
                                 * np.sin( X2[indX2] ), x = X2[indX2] )
                    G1[i,j] \
                      = simps( Data[j,indX2] * P1[indX2] \
                                 * np.sin( X2[indX2] ), x = X2[indX2] )
                    G2[i,j] \
                      = simps( Data[j,indX2] * P2[indX2] \
                                 * np.sin( X2[indX2] ), x = X2[indX2] )

                    G3[i,j] \
                      = simps( Data[j,indX2] * P3[indX2] \
                                 * np.sin( X2[indX2] ), x = X2[indX2] )

                    G4[i,j] \
                      = simps( Data[j,indX2] * P4[indX2] \
                                 * np.sin( X2[indX2] ), x = X2[indX2] )

                # --- Integrate over radial dimension ---

                if self.R0 < 0.0:

                    Power[i,0] \
                      = TwoPi \
                          * simps( G0[i,indX1]**2 * X1[indX1]**2, \
                                   x = X1[indX1] )

                    Power[i,1] \
                      = TwoPi \
                          * simps( G1[i,indX1]**2 * X1[indX1]**2, \
                                   x = X1[indX1] )

                    Power[i,2] \
                      = TwoPi \
                          * simps( G2[i,indX1]**2 * X1[indX1]**2, \
                                   x = X1[indX1] )

                    Power[i,3] \
                      = TwoPi \
                          * simps( G3[i,indX1]**2 * X1[indX1]**2, \
                                   x = X1[indX1] )

                    Power[i,4] \
                      = TwoPi \
                          * simps( G4[i,indX1]**2 * X1[indX1]**2, \
                                   x = X1[indX1] )

                else:

                    Power[i,0] = G0[i,indX1]**2
                    Power[i,1] = G1[i,indX1]**2
                    Power[i,2] = G2[i,indX1]**2
                    Power[i,3] = G3[i,indX1]**2
                    Power[i,4] = G4[i,indX1]**2

            Time0, RsAve, RsMin, RsMax = self.GetShockRadiusVsTime()

            Data = np.vstack( (Time,RsAve,RsMin,RsMax, \
                               Power[:,0],Power[:,1], \
                               Power[:,2],Power[:,3],Power[:,4]) )
            np.savetxt( self.FileName, Data )

        Time, RsAve, RsMin, RsMax, P0, P1, P2, P3, P4 \
          = np.loadtxt( self.FileName )

#        if( isfile( self.ShockRadiusVsTimeFileName ) ):
#            os.system( 'rm -f {:}'.format( self.ShockRadiusVsTimeFileName ) )

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

#        beta = self.CurveFit( tFit, np.log( P1[ind] ), InitialGuess )
#        F = np.exp( self.FittingFunction( tFit, beta ) )

        beta, pcov = curve_fit( FittingFunction2, tFit, np.log( P1[ind] ), \
                                p0 = InitialGuess, jac = Jacobian2 )
        F = np.exp( FittingFunction2( tFit, beta[0], beta[1], \
                                            beta[2], beta[3] ) )

        self.beta = np.copy( beta )
        self.perr = np.sqrt( np.diag( pcov ) )

        # Propagate error from frequency into period
        self.perr[1] = self.perr[1] / ( 2.0 * self.beta[1]**2 )
        self.perr[2] = TwoPi * self.perr[2] / self.beta[2]**2

        self.beta[0] = np.exp( self.beta[0] )
        self.beta[1] = 1.0 / ( 2.0 * self.beta[1] )
        self.beta[2] = TwoPi / self.beta[2]
        self.beta[3] = self.beta[3]

        if self.Verbose:
            print( '' )
            print( 'F1 = {:.3e}'.format   ( self.beta[0] ) )
            print( 'tr = {:.3e} ms'.format( self.beta[1] ) )
            print( 'ti = {:.3e} ms'.format( self.beta[2] ) )
            print( 'd  = {:.3e}'.format   ( self.beta[3] ) )

        return tFit+t0, F

    def PlotData \
          ( self, ax, t0, t1, \
            Time, RsAve, RsMin, RsMax, P0, P1, P2, P3, P4, tF, F, \
            M = -1.0, Rs = -1.0, GR = False ):

        ind = np.where( ( Time >= t0 ) & ( Time <= t1 ) )[0]

        t = np.copy( Time )

        Time  = Time [ind]
        RsAve = RsAve[ind]
        RsMin = RsMin[ind]
        RsMax = RsMax[ind]
        P0    = P0   [ind]
        P1    = P1   [ind]
        P2    = P2   [ind]
        P3    = P3   [ind]
        P4    = P4   [ind]

        if M > 0.0: ax.set_title( 'M{:}_Mdot0.3_Rs{:}'.format( M, Rs ) )

        c  = 'r'
        suffix = '(NR)'
        if GR:
            c = 'b'
            suffix = '(GR)'

        ax.plot( Time, P1, c + '-', label = 'P1 ' + suffix )

        if type( F ) == np.ndarray:

            ind = np.where( ( t >= tF[0] ) & ( t <= tF[-1] ) )[0]

            Time = t[ind]

            ax.plot( Time, F, c + '--', label = 'Fit ' + suffix )

#            txt = r'$\tau$ = ( {:.3e} $\pm$ {:.3e} ) ms'.format \
#                    ( self.beta[1], self.perr[1] )
#            txt += '\n'
#            txt += r'$T$ = ( {:.3e} $\pm$ {:.3e} ) ms'.format \
#                     ( self.beta[2], self.perr[2] )
#
#            ax.text( 0.2, 0.8, txt, transform = ax.transAxes )

        ax.set_yscale( 'log' )

        xlim = ( t0, t1 )
        ax.set_xlim( xlim )

        ax.set_xlabel( 'Time [ms]' )
        ax.set_ylabel( 'Power [cgs]' )

        ax.grid()

        return


if __name__ == "__main__":

    #Root = '/scratch/dunhamsj/ProductionRuns/'
    Root = '/lump/data/AccretionShockStudy/'

    Field = 'Entropy'
    t0    = 000.0
    t1    = 150.0
    fL    = 0.8
    fU    = 0.9
    R0    = -1.7e2
    suffix = ''

    M     = np.array( [ '1.4', '2.0', '2.8' ], str )
    Mdot  = '0.3'
    Rs    = np.array( [ '120', '150', '180' ], str )

    fig, axs = plt.subplots( 3, 3, figsize = (16,9) )

    for m in range( M.shape[0] ):
        for rs in range( Rs.shape[0] ):

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

            omega_r = 1.0 / ( 2.0 * tauR )
            omega_i = TwoPi / T_SASI

            InitialGuess = np.array( [ LogF0, omega_r, omega_i, delta ], \
                                     np.float64 )

            ID_NR = 'NR2D_M{:}_Mdot{:}_Rs{:}'.format( M[m], Mdot, Rs[rs] )
            P_NR = PowersInLegendreModes( Root, ID_NR, Field, \
                                          Rs = np.float64( Rs[rs] ), \
                                          fL = fL, fU = fU, R0 = R0, \
                                          EntropyThreshold = 4.0e14, \
                                          Verbose = False )
            Time, RsAve, RsMin, RsMax, P0, P1, P2, P3, P4 \
              = P_NR.ComputePowerInLegendreModes()
            tFit, F = P_NR.FitPowerInLegendreModes \
                         ( Time, tF0, tF1, P1, InitialGuess = InitialGuess )
            P_NR.PlotData( axs[m,rs], t0, t1, Time, RsAve, RsMin, RsMax, \
                           P0, P1, P2, P3, P4, tFit, F )
            del ID_NR, P_NR, Time, RsAve, RsMin, RsMax, \
                P0, P1, P2, P3, P4, tFit, F

            ID_GR = 'GR2D_M{:}_Mdot{:}_Rs{:}'.format( M[m], Mdot, Rs[rs] )
            P_GR = PowersInLegendreModes( Root, ID_GR, Field, \
                                          Rs = np.float64( Rs[rs] ), \
                                          fL = fL, fU = fU, R0 = R0, \
                                          EntropyThreshold = 4.0e14, \
                                          Verbose = False )
            Time, RsAve, RsMin, RsMax, P0, P1, P2, P3, P4 \
              = P_GR.ComputePowerInLegendreModes()
            tFit, F = P_GR.FitPowerInLegendreModes \
                        ( Time, tF0, tF1, P1, InitialGuess = InitialGuess )
            P_GR.PlotData( axs[m,rs], t0, t1, Time, RsAve, RsMin, RsMax, \
                           P0, P1, P2, P3, P4, tFit, F, \
                           np.float64( M[m] ), np.float64( Rs[rs] ), GR = True )
            del ID_GR, P_GR, Time, RsAve, RsMin, RsMax, \
                P0, P1, P2, P3, P4, tFit, F

    axs[0,0].legend()
    plt.subplots_adjust( wspace = 0.3, hspace = 0.4 )
    plt.savefig( \
    '/home/kkadoogan/fig.PowersInLegendreModes_MultiRun_{:}.png'.format \
    ( Field ), dpi = 300 )

    import os
    os.system( 'rm -rf __pycache__ ' )
