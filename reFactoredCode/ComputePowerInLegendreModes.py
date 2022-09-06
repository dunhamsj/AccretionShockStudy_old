#!/usr/bin/env python3

from scipy.integrate import simps
import numpy as np

from UtilitiesModule import GetFileArray, Overwrite, ComputeAngleAverage
from PowerSpectrumUtilitiesModule import ReadFields

def ComputePowerInLegendreModes \
      ( plotFileDirectory, plotFileBaseName, dataFileName, \
        field, fL, fU, Rs, verbose = False ):

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

    OW = Overwrite( dataFileName )

    if( not OW ): return

    plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )

    nSS = plotFileArray.shape[0]

    H = np.empty( (nSS,5), np.float64 )

    time = np.empty( nSS )

    nLeg = 5

    Rs *= 1.0e5

    for iSS in range( nSS ):

        if( verbose ):
            if( ( iSS + 1 ) % 10 == 0 ):
              print( '  File {:}/{:}'.format( iSS + 1, nSS ) )

        plotFile = plotFileDirectory + plotFileArray[iSS]

        time[iSS], data, X1, X2, X3, dX1, dX2, dX3, nX \
          = ReadFields( plotFile, field )

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

        G = np.empty( (nLeg,nX[0]), np.float64 )

        A = np.zeros( (nX[0],nX[1],nX[2]), np.float64 )

        if( field == 'DivV2' ):

            # --- Sheck et al., (2008), A&A, 477, 931 ---

            PF_V2 = data[1]

            indX1 = np.where( ( X1 > fL * Rs ) & ( X1 < fU * Rs ) )[0]
            indX2 = np.linspace( 1, nX[1]-2, nX[1]-2, dtype = np.int64 )
            indX3 = np.linspace( 0, nX[2]  , nX[2]  , dtype = np.int64 )

            for i in indX1:
                for j in indX2:
                    for k in indX3:
                        A[i,j,k] \
                          = 1.0 / ( 2.0 * dX2[j] * np.sin( X2[j] ) ) \
                              * (   np.sin( X2[j+1] ) * PF_V2[i,j+1,k] \
                                  - np.sin( X2[j-1] ) * PF_V2[i,j-1,k] )

        else:

            print( 'Invalid choice of field: {:}'.format( Field ) )
            print( 'Valid choices' )
            print( '-------------' )
            print( '  DivV2' )
            exit( 'Exiting...' )

        Psi    = data[0]
        Psi_AA = np.empty( nX[0], np.float64 )

        for i in indX1:

          Psi_AA[i] \
            = ComputeAngleAverage \
                ( Psi[i,indX2,indX3], X2[indX2], dX2[indX2], dX3[indX3] )

          # Subtract off angle-average
          A[i,indX2,indX3] \
            -= ComputeAngleAverage \
                 ( A[i,indX2,indX3], X2[indX2], dX2[indX2], dX3[indX3] )

        #Psi_AA[indX1] = 1.0

        for p in range( nLeg ):

            # --- Compute p-th expansion coefficient ---

            for i in indX1:

                G[p,i] \
                  = simps( A[i,indX2,indX3] * P[p,indX2] \
                             * np.sin( X2[indX2] ), x = X2[indX2] )

            # --- Integrate over radial dimension ---

            H[iSS,p] \
              = 2.0 * np.pi \
                  * simps( G[p,indX1]**2 * Psi_AA[indX1]**6 * X1[indX1]**2, \
                           x = X1[indX1] )

    # END for iSS in range( nSS )

    Data = np.vstack( (time,H[:,0],H[:,1],H[:,2],H[:,3],H[:,4]) )
    np.savetxt( dataFileName, Data )

    return
#END ComputePowerInLegendreModes

if __name__ == '__main__':

    ID = 'GR2D_M2.8_Mdot0.3_Rs120'
    plotFileDirectory \
      = '/lump/data/accretionShockStudy/{:}/'.format( ID )
    plotFileBaseName = '{:}.plt_'.format( ID )
    dataFileName = '.LegendrePowerSpectrum_{:}'.format( ID )
    ComputePowerInLegendreModes \
      ( plotFileDirectory, plotFileBaseName, dataFileName, \
        'DivV2', 0.8, 0.9, 120.0, verbose = True )
