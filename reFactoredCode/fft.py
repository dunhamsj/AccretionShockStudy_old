#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from FitPowerToModel import FitPowerToModel, FittingFunction

ID = '2D_M2.8_Mdot0.3_Rs6.00e1_RPNS2.00e1'
#ID = '2D_M2.8_Mdot0.3_Rs120'

ID_GR = 'GR' + ID
ID_NR = 'NR' + ID
timeGR, dataGR = np.loadtxt( 'LatFlux_{:}.dat'.format( ID_GR ) )
timeNR, dataNR = np.loadtxt( 'LatFlux_{:}.dat'.format( ID_NR ) )

dataFileName = 'dataFileName.dat'

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
FitPowerToModel( t0, t1, timeGR, dataGR**2, InitialGuess, dataFileName )

t0_GR, t1_GR, LogF_GR, omegaR_GR, omegaI_GR, delta_GR, \
dummy0, omegaR_err_GR, omegaI_err_GR, dummy1 \
  = np.loadtxt( dataFileName )

ind_GR = np.where( ( timeGR >= t0_GR ) & ( timeGR <= t1_GR ) )[0]

tF_GR = timeGR[ind_GR]

F_GR = FittingFunction \
         ( tF_GR - tF_GR[0], LogF_GR, omegaR_GR, omegaI_GR, delta_GR )

fig, ax = plt.subplots( 1, 1 )
fig.suptitle( ID_GR )

ax.semilogy( timeGR, dataGR**2 )
ax.semilogy( tF_GR, np.exp( F_GR ) )
ax.set_xlabel( 'Time [ms]' )

#plt.show()
plt.savefig( '/home/kkadoogan/power60.png', dpi = 300 )
plt.close()

fig, axs = plt.subplots( 2, 1 )
fig.suptitle( ID_GR )

axs[0].plot( timeGR, dataGR )
axs[0].set_xlabel( 'Time [ms]' )
axs[0].set_ylabel( 'data' )

N = timeGR.shape[0]
A = np.abs( np.fft.fft( dataGR )[1:N//2] )**2
f = np.fft.fftfreq( N )[1:N//2]
axs[1].plot( 1.0 / f, A.real )
axs[1].set_xlabel( 'Period [ms]' )
axs[1].set_ylabel( 'FFT Amplitude' )

plt.subplots_adjust( hspace = 0.3 )
plt.show()
#plt.savefig( '/home/kkadoogan/fft.png', dpi = 300 )
plt.close()
