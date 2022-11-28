#!/usr/bin/env python3

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

ID = 'GR2D_M2.8_Mdot0.3_Rs180'

timeA, dataA = np.loadtxt( 'LatFlux_{:}.dat'.format( ID ) )

indmax = timeA.shape[0]

indA = np.linspace( 0, indmax-1, indmax, dtype = np.int64 )

fig, ax = plt.subplots( 1, 1 )
ax.set_title( r'\texttt{{{:}}}'.format( ID ) )

t = np.linspace( 100.0, 300.0, 5 )

colors = plt.cm.inferno( np.linspace( 0.25, 0.75, t.shape[0] ) )
alpha = np.linspace( 0.2, 1.0, t.shape[0] )

for i in range( t.shape[0] ):

    ind = np.where( timeA < t[i] )[0]

    time = np.copy( timeA[ind] )
    data = np.copy( dataA[ind] )

    # Compute and plot FFT

    N = time.shape[0]

    if( N % 2 == 0 ):
        y = np.abs( sp.fft.fft    ( data )[1:N//2] )**2
        x =         sp.fft.fftfreq( N, 1 )[1:N//2]
    else:
        y = np.abs( sp.fft.fft    ( data )[1:(N-1)//2] )**2
        x =         sp.fft.fftfreq( N, 1 )[1:(N-1)//2]

    x = 1.0 / x

    x = x[::-1]
    y = y[::-1]

    ax.plot( x, np.abs( y / y.max() ), \
            #             color = colors[i], \
             'k-', \
             alpha = alpha[i], \
             label = r'$T_{{\mathrm{{SASI}}}}$: {:.2f} ms' \
                     .format( x[y.argmax()] ) )

ax.legend()
ax.set_xlabel( 'Time [ms]' )
ax.set_ylabel( 'Normalized FFT Amplitude' )

ax.grid()

plt.savefig( 'fig.{:}_FFT_Normed.png'.format( ID ), dpi = 300 )
#plt.show()
plt.close()
