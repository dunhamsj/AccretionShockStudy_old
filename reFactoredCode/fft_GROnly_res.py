#!/usr/bin/env python3

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

ID = 'GR2D_M2.8_Mdot0.3_Rs6.00e1_RPNS2.00e1'

timeA, dataA = np.loadtxt( 'LatFlux_{:}_highRes.dat'.format( ID ) )

indmax = timeA.shape[0]

indA = np.linspace( 0, indmax-1, indmax, dtype = np.int64 )

fig, ax = plt.subplots( 1, 1 )
ax.set_title( r'\texttt{{{:}}}'.format( ID ) )

nE = np.array( [ 1, 2, 4, 8 ], dtype = np.int64 )
alpha = np.linspace( 0.2, 1.0, nE.shape[0] )

for i in range( nE.shape[0] ):

    time = np.copy( timeA[::nE[i]] )
    data = np.copy( dataA[::nE[i]] )

    # Compute and plot FFT

    N = time.shape[0]

    T = time[-1] / N

    if( N % 2 == 0 ):
        y = np.abs( sp.fft.fft    ( data )[1:N//2] )**2
        x =         sp.fft.fftfreq( N, T )[1:N//2]
    else:
        y = np.abs( sp.fft.fft    ( data )[1:(N-1)//2] )**2
        x =         sp.fft.fftfreq( N, T )[1:(N-1)//2]

    x = 1.0 / x

    x = x[::-1]
    y = y[::-1]

    ax.plot( x, np.abs( y / y.max() ), \
             'k-', \
             alpha = alpha[i], \
             label = r'nE: {:}, $T_{{\mathrm{{SASI}}}}$: {:.2f} ms' \
                     .format( nE[i], x[y.argmax()] ) )

ax.legend()
ax.set_xlabel( 'Time [ms]' )
ax.set_ylabel( 'Normalized FFT Amplitude' )

ax.grid()

plt.savefig( 'fig.{:}_FFT_Normed_skipEvery.png'.format( ID ), dpi = 300 )
#plt.show()
plt.close()
