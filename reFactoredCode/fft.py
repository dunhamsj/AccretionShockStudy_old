#!/usr/bin/env python3

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from computeTimeScales import ComputeTimeScales

M    = '2.8'
Rs   = '7.00e1'
Rpns = '020'
ID = '2D_M{:}_Rpns{:}_Rs{:}'.format( M, Rpns, Rs )

ID_GR = 'GR' + ID
ID_NR = 'NR' + ID
timeGR, dataGR = np.loadtxt( 'LatFlux_{:}.dat'.format( ID_GR ) )
timeNR, dataNR = np.loadtxt( 'LatFlux_{:}.dat'.format( ID_NR ) )

indmax = min( timeGR.shape[0], timeNR.shape[0] )
#ind = np.linspace( 0, indmax-1, indmax, dtype = np.int64 )
ind = -1#461 # Rs = 8.75e1
#ind = 195  # Rs = 6.00e1
#ind = 777  # Rs = 1.20e2
#ind = 1330 # Rs = 1.50e2
#ind = 1839 # Rs = 1.75e2

timeGR = np.copy( timeGR[0:ind] )
timeNR = np.copy( timeNR[0:ind] )
dataGR = np.copy( dataGR[0:ind] )
dataNR = np.copy( dataNR[0:ind] )

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

rInner = np.float64( Rpns )
rOuter = np.float64( Rs )

plotFileDirectory = '/lump/data/accretionShockStudy/newData/2D/'

plotFileDirectory_NR = plotFileDirectory + ID_NR + '/'
plotFileDirectory_NR += ID_NR + '.plt00000000/'
tauAd, tauAc = ComputeTimeScales( plotFileDirectory_NR, rInner, rOuter, 'NR' )
T_SASI_NR = tauAd + tauAc

plotFileDirectory_GR = plotFileDirectory + ID_GR + '/'
plotFileDirectory_GR += ID_GR + '.plt00000000/'
tauAd, tauAc = ComputeTimeScales( plotFileDirectory_GR, rInner, rOuter, 'GR' )
T_SASI_GR = tauAd + tauAc

print()
print( 'T_SASI' )
print( 'Muller (NR): {:.3e} ms'.format( T_SASI_NR ) )
print( 'FFT    (NR): {:.3e} ms'.format( xNR[yNR.argmax()] ) )
print( '|deltaT/T| (NR): {:.3e}' \
       .format( abs( ( T_SASI_NR - xNR[yNR.argmax()] ) / T_SASI_NR ) ) )
print( 'Muller (GR): {:.3e} ms'.format( T_SASI_GR ) )
print( 'FFT    (GR): {:.3e} ms'.format( xGR[yGR.argmax()] ) )
print( '|deltaT/T| (GR): {:.3e}' \
       .format( abs( ( T_SASI_GR - xNR[yGR.argmax()] ) / T_SASI_GR ) ) )
print( 'Mulller (GR/NR): {:.3e}'.format( T_SASI_GR / T_SASI_NR ) )
print( 'FFT     (GR/NR): {:.3e}'.format( xGR[yGR.argmax()] / xNR[yNR.argmax()] ) )

# colorblind-friendly palette: https://gist.github.com/thriveth/8560036
color = ['#377eb8', '#ff7f00', '#4daf4a', \
         '#f781bf', '#a65628', '#984ea3', \
         '#999999', '#e41a1c', '#dede00']

fig, axs = plt.subplots( 2, 1 )
fig.suptitle( r'\texttt{{{:}}}'.format( ID ) )

axs[0].plot( timeNR, dataNR, color = color[0], label = 'NR' )
axs[0].plot( timeGR, dataGR, color = color[1], label = 'GR' )

axs[1].text( 0.5, 0.6, 'GR: {:.3e} ms'.format( xGR[yGR.argmax()] ), \
             transform = axs[1].transAxes )
axs[1].text( 0.5, 0.5, 'NR: {:.3e} ms'.format( xNR[yNR.argmax()] ), \
             transform = axs[1].transAxes )

axs[1].plot( xNR, np.abs( yNR / yNR.max() ), '-', color = color[0] )
axs[1].plot( xGR, np.abs( yGR / yGR.max() ), '-', color = color[1] )

#axs[0].set_xlim( 0.0, 20)#1.0e2 )
#axs[1].set_xlim( 0.0, 20)#1.0e2 )
axs[0].set_xlabel( 'Time [ms]' )
axs[1].set_xlabel( 'Time [ms]' )
axs[0].set_ylabel( 'Data' )
axs[1].set_ylabel( 'Normalized FFT Amplitude' )

axs[0].legend()

axs[0].grid()
axs[1].grid()

plt.subplots_adjust( hspace = 0.3 )
#plt.savefig( '/home/kkadoogan/fig.{:}_FFT_Normed.png'.format( ID ), dpi = 300 )
plt.show()
plt.close()
