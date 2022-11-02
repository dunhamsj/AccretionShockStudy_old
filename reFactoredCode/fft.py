#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from FitPowerToModel import FitPowerToModel, FittingFunction

#ID = '2D_M2.8_Mdot0.3_Rs6.00e1_RPNS2.00e1'
ID = '2D_M2.8_Mdot0.3_Rs120'

ID_GR = 'GR' + ID
ID_NR = 'NR' + ID
timeGR, dataGR = np.loadtxt( 'LatFlux_{:}.dat'.format( ID_GR ) )
timeNR, dataNR = np.loadtxt( 'LatFlux_{:}.dat'.format( ID_NR ) )

dataFileNameGR = '{:}_FFT.dat'.format( ID_GR )
dataFileNameNR = '{:}_FFT.dat'.format( ID_NR )

t0 = 3.0
t1 = 7.5e1

LogF  = np.log( 1.0e14 )
tauR  = 200.0
delta = 0.0
T_SASI = 9.1

omegaR = 2.0 * np.pi / tauR
omegaI = 2.0 * np.pi / T_SASI

InitialGuess \
  = np.array( [ LogF, omegaR, omegaI, delta ], np.float64 )
FitPowerToModel( t0, t1, timeGR, dataGR**2, InitialGuess, dataFileNameGR )
FitPowerToModel( t0, t1, timeNR, dataNR**2, InitialGuess, dataFileNameNR )

t0_GR, t1_GR, LogF_GR, omegaR_GR, omegaI_GR, delta_GR, \
dummy0, omegaR_err_GR, omegaI_err_GR, dummy1 \
  = np.loadtxt( dataFileNameGR )

t0_NR, t1_NR, LogF_NR, omegaR_NR, omegaI_NR, delta_NR, \
dummy0, omegaR_err_NR, omegaI_err_NR, dummy1 \
  = np.loadtxt( dataFileNameNR )

ind_GR = np.where( ( timeGR >= t0_GR ) & ( timeGR <= t1_GR ) )[0]
ind_NR = np.where( ( timeNR >= t0_NR ) & ( timeNR <= t1_NR ) )[0]

tF_GR = timeGR[ind_GR]
tF_NR = timeNR[ind_NR]

F_GR = FittingFunction \
         ( tF_GR - tF_GR[0], LogF_GR, omegaR_GR, omegaI_GR, delta_GR )
F_NR = FittingFunction \
         ( tF_NR - tF_NR[0], LogF_NR, omegaR_NR, omegaI_NR, delta_NR )

fig, ax = plt.subplots( 1, 1 )
fig.suptitle( ID )

ax.semilogy( timeGR, dataGR**2 )
ax.semilogy( timeNR, dataNR**2 )
ax.semilogy( tF_GR, np.exp( F_GR ) )
ax.semilogy( tF_NR, np.exp( F_NR ) )
ax.set_xlabel( 'Time [ms]' )

plt.show()
#plt.savefig( '/home/kkadoogan/power60.png', dpi = 300 )
plt.close()

fig, axs = plt.subplots( 2, 1 )
fig.suptitle( ID )

axs[0].plot( timeGR, dataGR )
axs[0].plot( timeNR, dataNR )
axs[0].set_xlabel( 'Time [ms]' )
axs[0].set_ylabel( 'data' )

NGR = timeGR.shape[0]
AGR = np.abs( np.fft.fft( dataGR )[1:NGR//2] )**2
fGR = np.fft.fftfreq( NGR )[1:NGR//2]
axs[1].plot( 1.0 / fGR, AGR.real )
NNR = timeNR.shape[0]
ANR = np.abs( np.fft.fft( dataNR )[1:NNR//2] )**2
fNR = np.fft.fftfreq( NNR )[1:NNR//2]
axs[1].plot( 1.0 / fNR, ANR.real )

axs[1].set_xlabel( 'Period [ms]' )
axs[1].set_ylabel( 'FFT Amplitude' )
axs[1].set_yscale( 'log' )

plt.subplots_adjust( hspace = 0.3 )
plt.show()
#plt.savefig( '/home/kkadoogan/fft.png', dpi = 300 )
plt.close()
