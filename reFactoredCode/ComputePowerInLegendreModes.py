#!/usr/bin/env python3

from scipy.integrate import trapezoid
import numpy as np

from UtilitiesModule import GetFileArray, Overwrite
from PowerSpectrumUtilitiesModule import ReadFields

def ComputePowerInLegendreModes \
      ( plotFileDirectory, plotFileBaseName, dataFileName, \
        field, fL, fU, Rs, fc = False, ow = False, verbose = False ):

    if( verbose ):
        print()
        print( '  Calling ComputePowerInLegendreModes...' )
        print( '  --------------------------------------' )
        print()
        print( '{:>21} : {:}'.format( 'plotFileDirectory', plotFileDirectory ) )
        print( '{:>21} : {:}'.format( 'plotFileBaseName', plotFileBaseName ) )
        print( '{:>21} : {:}'.format( 'dataFileName', dataFileName ) )
        print( '{:>21} : {:}'.format( 'field', field ) )
        print( '{:>21} : {:}'.format( 'fL', fL ) )
        print( '{:>21} : {:}'.format( 'fU', fU ) )
        print( '{:>21} : {:} km'.format( 'Rs', Rs ) )
        print()

    OW = Overwrite( dataFileName, ForceChoice = fc, OW = ow )

    if( not OW ): return

    plotFileArray \
      = GetFileArray( plotFileDirectory, plotFileBaseName )
    plotFileArray = np.copy( plotFileArray[:-1] ) # remove 99999999 file

    nSS = plotFileArray.shape[0]

    H = np.zeros( (nSS,5), np.float64 )

    time = np.empty( nSS )

    nLeg = 5

    Rs *= 1.0e5

    for iSS in range( nSS ):

        plotFile = plotFileDirectory + plotFileArray[iSS]

        time[iSS], data, X1, X2, X3, dX1, dX2, dX3, nX \
          = ReadFields( plotFile, field )

        if( verbose ):
            if( ( iSS + 1 ) % 10 == 0 ):
                print( '  File {:}/{:}, t: {:.6f} ms' \
                       .format( iSS + 1, nSS, time[iSS] ) )

        X1  *= 1.0e5
        dX1 *= 1.0e5

        x = np.cos( X2 )

        # Legendre polynomials normalized s.t. integral( P_m, P_n ) = delta_mn

        P = np.empty( (nLeg,nX[1]), np.float64 )

        P[0] = np.sqrt( 1.0 / 2.0 ) \
                 * np.ones( nX[1] )
        P[1] = np.sqrt( 3.0 / 2.0 ) \
                 * x
        P[2] = np.sqrt( 5.0 / 2.0 ) \
                 * ( 3.0 * x**2 - 1.0 ) / 2.0
        P[3] = np.sqrt( 7.0 / 2.0 ) \
                 * 1.0 / 2.0 * ( 5.0 * x**3 - 3.0 * x )
        P[4] = np.sqrt( 9.0 / 2.0 ) \
                 * 1.0 / 8.0 * ( 35.0 * x**4 - 30.0 * x**2 + 3.0 )

        G = np.zeros( (nLeg,nX[0]), np.float64 )

        A = np.zeros( (nX[0],nX[1],nX[2]), np.float64 )

        if( field == 'DivV2' ):

            # --- Sheck et al., (2008), A&A, 477, 931 ---

            PF_V2 = data[1]

            indX1 = np.where( ( X1 > fL * Rs ) & ( X1 < fU * Rs ) )[0]
            indX2 = np.linspace( 0, nX[1]-1, nX[1], dtype = np.int64 )
            indX3 = np.linspace( 0, nX[2]  , nX[2], dtype = np.int64 )

            for i in indX1:
                for j in indX2:
                    for k in indX3:

                        # Reflecting boundary conditions in theta
                        if( j == 0 ):
                            X2m = X2[j]
                            X2p = X2[j+1]
                            V2m = -PF_V2[i,j  ,k]
                            V2p = +PF_V2[i,j+1,k]
                        elif( j == nX[1]-1 ):
                            X2m = X2[j-1]
                            X2p = X2[j]
                            V2m = +PF_V2[i,j-1,k]
                            V2p = -PF_V2[i,j  ,k]
                        else:
                            X2m = X2[j-1]
                            X2p = X2[j+1]
                            V2m = PF_V2[i,j-1,k]
                            V2p = PF_V2[i,j+1,k]

                        A[i,j,k] \
                          = 1.0 / ( 2.0 * dX2[j] * np.sin( X2[j] ) ) \
                              * (   np.sin( X2p ) * V2p \
                                  - np.sin( X2m ) * V2m )

        elif( field == 'PolytropicConstant' ):

            indX1 = np.where( ( X1 > fL * Rs ) & ( X1 < fU * Rs ) )[0]
            indX2 = np.linspace( 0, nX[1]-1, nX[1], dtype = np.int64 )
            indX3 = np.linspace( 0, nX[2]  , nX[2], dtype = np.int64 )

            A = data[1]

        else:

            print( 'Invalid choice of field: {:}'.format( field ) )
            print( 'Valid choices' )
            print( '-------------' )
            print( '  DivV2' )
            print( '  PolytropicConstant' )
            exit( 'Exiting...' )

        Psi    = data[0]
        Psi_AA = np.empty( nX[0], np.float64 )

#        for i in indX1:
#
#          Psi_AA[i] \
#            = ComputeAngleAverage \
#                ( Psi[i,indX2,indX3], X2[indX2], dX2[indX2], dX3[indX3] )
#
#          # Subtract off angle-average
#          A[i,indX2,indX3] \
#            -= ComputeAngleAverage \
#                 ( A[i,indX2,indX3], X2[indX2], dX2[indX2], dX3[indX3] )

        #Psi_AA[indX1] = 1.0

        for ell in range( nLeg ):

            # --- Compute ell-th expansion coefficient ---

            for i in indX1:

                # Assume data is 2D
                for k in indX3:

                    G[ell,i] \
                      += trapezoid \
                           ( A[i,indX2,k] * P[ell,indX2] \
                               * np.sin( X2[indX2] ), \
                             x = X2[indX2], dx = dX2[0] )

            # Assume Psi is spherically symmetric
            H[iSS,ell] \
              = 4.0 * np.pi \
                  * trapezoid \
                      ( G[ell,indX1]**2 * Psi[indX1,0,0]**6 * X1[indX1]**2, \
                        x = X1[indX1], dx = dX1[0] )

    # END for iSS in range( nSS )

    Data = np.vstack( (time,H[:,0],H[:,1],H[:,2],H[:,3],H[:,4]) )
    np.savetxt( dataFileName, Data )

    return
#END ComputePowerInLegendreModes

if __name__ == '__main__':

    ID = 'NR2D_M2.8_Mdot0.3_Rs120_nX768x064'
    plotFileDirectory \
      = '/lump/data/accretionShockStudy/angularResolution/{:}/'.format( ID )
    plotFileBaseName = '{:}.plt_'.format( ID )
    dataFileName = '.{:}_LegendrePowerSpectrum.dat'.format( ID )
    ComputePowerInLegendreModes \
      ( plotFileDirectory, plotFileBaseName, dataFileName, \
        'DivV2', 0.8, 0.9, 120.0, verbose = True )
