#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from UtilitiesModule import GetData

Data \
  = np.loadtxt( '.GR2D_M1.4_Mdot0.3_Rs180_PA0.004_nX320x064_smooth_DivV2_110.00km_PowersInLegendreModes.dat' )

Time = Data[0]
P1   = Data[2]

ind = np.where( Time > 20.0 )[0]

Time = Time[ind]
P1   = np.log10( P1[ind] )

P1k = np.fft.fft( P1 )
P1t = np.fft.ifft( P1k )

fig, axs = plt.subplots( 2, 1 )

axs[0].plot( Time, P1 , '.', label = 'Orig' )
axs[0].plot( Time, P1t, 'x', label = 'invFFT' )
axs[0].legend()
axs[0].set_xlabel( 'Time [ms]' )
axs[0].set_ylabel( 'log of Power in P1 [cgs]' )

P1k = np.abs( P1k )
freq = np.fft.fftfreq( Time.shape[-1] )

ind = np.where( freq > 0.0 )[0]

freq = freq[ind]
P1k  = P1k[ind]

axs[1].plot( freq, P1k )
axs[1].set_xlabel( 'Frequency [Hz]' )
axs[1].set_ylabel( 'Amplitude' )
axs[1].set_yscale( 'log' )
axs[1].axvline( 2.0 * np.pi / 50.0, c='k' )
axs[1].text( 0.13, 1.0e2, '50 ms' )

plt.subplots_adjust( hspace = 0.5 )
#plt.show()
plt.savefig( 'fig.FFT.png', dpi = 300 )
