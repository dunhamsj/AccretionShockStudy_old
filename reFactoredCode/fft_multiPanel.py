#!/usr/bin/env python3

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from computeTimeScales import ComputeTimeScales

Rs   = [ '6.00e1', '7.50e1', '9.00e1' ]
RPNS = '2.00e1'

cutoff = 0.03
fig, axs = plt.subplots( 3, 1 )
fig.suptitle( r'$M=2.8\ \mathrm{M}_{\odot}$' + '\n' \
              + r'$\left(R_{{s}}-R_{{s,0}}\right)/R_{{s,0}} < {:}$' \
                .format( cutoff ) )
fig.suptitle( r'$M=2.8\ \mathrm{M}_{\odot}$' )

for rs in range( len( Rs ) ):

    ID = '2D_M2.8_Mdot0.3_Rs{:}_RPNS{:}'.format( Rs[rs], RPNS )

    ID_GR = 'GR' + ID
    ID_NR = 'NR' + ID
    timeGR, dataGR = np.loadtxt( 'LatFlux_{:}.dat'.format( ID_GR ) )
    timeNR, dataNR = np.loadtxt( 'LatFlux_{:}.dat'.format( ID_NR ) )

    RsFileNameGR = '{:}_ShockRadiusVsTime.dat'.format( ID_GR )
    TimeGR, RsAveGR, RsMinGR, RsMaxGR = np.loadtxt( RsFileNameGR )
    RsFileNameNR = '{:}_ShockRadiusVsTime.dat'.format( ID_NR )
    TimeNR, RsAveNR, RsMinNR, RsMaxNR = np.loadtxt( RsFileNameNR )

    indmaxGR = timeGR.shape[0]
    indmaxNR = timeNR.shape[0]

    indGR = np.where( ( RsAveGR - RsAveGR[0] ) / RsAveGR[0] > cutoff )[0]
    indNR = np.where( ( RsAveNR - RsAveNR[0] ) / RsAveNR[0] > cutoff )[0]

    if indGR.shape[0] == 0:
        indGR = np.linspace( 0, indmaxGR-1, indmaxGR, dtype = np.int64 )

    if indNR.shape[0] == 0:
        indNR = np.linspace( 0, indmaxNR-1, indmaxNR, dtype = np.int64 )

    #indGR = np.linspace( 0, indmaxGR-1, indmaxGR, dtype = np.int64 )
    #indNR = np.linspace( 0, indmaxNR-1, indmaxNR, dtype = np.int64 )

    indGR = np.where( timeGR < 2.0e2 )[0]
    indNR = np.where( timeNR < 2.0e2 )[0]
    if Rs[rs] == '6.00e1':
        indNR = np.where( timeNR < 5.0e1 )[0]
    if Rs[rs] == '7.50e1':
        indNR = np.where( timeNR < 9.0e1 )[0]

    timeGR = np.copy( timeGR[indGR] )
    timeNR = np.copy( timeNR[indNR] )
    dataGR = np.copy( dataGR[indGR] )
    dataNR = np.copy( dataNR[indNR] )

    dataFileNameGR = '{:}_FFT.dat'.format( ID_GR )
    dataFileNameNR = '{:}_FFT.dat'.format( ID_NR )

    # Compute and plot FFT

    NGR = timeGR.shape[0]
    NNR = timeNR.shape[0]
    TGR = ( timeGR[-1] - timeGR[0] )/ np.float64( NGR ) # milliseconds
    TNR = ( timeNR[-1] - timeNR[0] )/ np.float64( NNR ) # milliseconds

    if( NGR % 2 == 0 ):
        yGR = np.abs( sp.fft.fft    ( dataGR )[1:NGR//2] )**2
        xGR =         sp.fft.fftfreq( NGR, TGR )[1:NGR//2]
    else:
        yGR = np.abs( sp.fft.fft    ( dataGR )[1:(NGR-1)//2] )**2
        xGR =         sp.fft.fftfreq( NGR, TGR )[1:(NGR-1)//2]

    if( NNR % 2 == 0 ):
        yNR = np.abs( sp.fft.fft    ( dataNR )[1:NNR//2] )**2
        xNR =         sp.fft.fftfreq( NNR, TNR )[1:NNR//2]
    else:
        yNR = np.abs( sp.fft.fft    ( dataNR )[1:(NNR-1)//2] )**2
        xNR =         sp.fft.fftfreq( NNR, TNR )[1:(NNR-1)//2]

    xGR = 1.0 / xGR
    xNR = 1.0 / xNR

    xGR = xGR[::-1]
    yGR = yGR[::-1]
    xNR = xNR[::-1]
    yNR = yNR[::-1]

    rInner = np.float64( RPNS )
    rOuter = np.float64( Rs[rs] )

    plotFileDirectory = '/lump/data/accretionShockStudy/'

    plotFileDirectory_NR = plotFileDirectory + ID_NR + '/'
    plotFileDirectory_NR += ID_NR + '.plt00000000/'
    tauAd, tauAc \
      = ComputeTimeScales( plotFileDirectory_NR, rInner, rOuter, 'NR' )
    T_SASI_NR = tauAd + tauAc

    plotFileDirectory_GR = plotFileDirectory + ID_GR + '/'
    plotFileDirectory_GR += ID_GR + '.plt00000000/'
    tauAd, tauAc \
      = ComputeTimeScales( plotFileDirectory_GR, rInner, rOuter, 'GR' )
    T_SASI_GR = tauAd + tauAc

    #print()
    #print( 'Rs = {:} km'.format( Rs[rs] ) )
    #print( 'Muller (NR): {:.3e} ms'.format( T_SASI_NR ) )
    #print( 'FFT    (NR): {:.3e} ms'.format( xNR[yNR.argmax()] ) )
    #print( '|deltaT/T| (NR): {:.3e}' \
    #       .format( abs( ( T_SASI_NR - xNR[yNR.argmax()] ) / T_SASI_NR ) ) )
    #print( 'Muller (GR): {:.3e} ms'.format( T_SASI_GR ) )
    #print( 'FFT    (GR): {:.3e} ms'.format( xGR[yGR.argmax()] ) )
    #print( '|deltaT/T| (GR): {:.3e}' \
    #       .format( abs( ( T_SASI_GR - xGR[yGR.argmax()] ) / T_SASI_GR ) ) )
    #print( 'Mulller (GR/NR): {:.3e}'.format( T_SASI_GR / T_SASI_NR ) )
    #print( 'FFT     (GR/NR): {:.3e}' \
    #       .format( xGR[yGR.argmax()] / xNR[yNR.argmax()] ) )

    axs[rs].text( 0.5, 0.8, \
                  r'$R_{{\mathrm{{s}}}}={:}\ \mathrm{{km}}$'.format( Rs[rs] ), \
                  transform = axs[rs].transAxes )
    axs[rs].text( 0.5, 0.6, 'GR: {:.3e} ms'.format( xGR[yGR.argmax()] ), \
                  transform = axs[rs].transAxes )
    axs[rs].text( 0.5, 0.5, 'NR: {:.3e} ms'.format( xNR[yNR.argmax()] ), \
                  transform = axs[rs].transAxes )

    axs[rs].plot( xGR, np.abs( yGR / yGR.max() ), '-', label = 'GR' )
    axs[rs].plot( xNR, np.abs( yNR / yNR.max() ), '-', label = 'NR' )

    axs[rs].set_xlabel( r'$T_{\mathrm{SASI}}\ \left[\mathrm{ms}\right]$' )

    axs[rs].grid()

axs[0].legend()
fig.supylabel( 'Normalized FFT Amplitude' )
plt.subplots_adjust( hspace = 0.0 )
plt.show()
#plt.savefig( '/home/kkadoogan/fig.FFT_LateStage.png', dpi = 300 )
#plt.savefig( '/home/kkadoogan/fig.FFT_LateStage_cutoff_{:}.png'.format(cutoff), dpi = 300 )
plt.close()
