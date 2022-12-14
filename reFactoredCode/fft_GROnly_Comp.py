#!/usr/bin/env python3

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

ID = 'GR2D_M2.8_Mdot0.3_Rs6.00e1_RPNS2.00e1'

timeA, dataA = np.loadtxt( 'LatFlux_{:}.dat'.format( ID ) )
timeB, dataB = np.loadtxt( 'LatFlux_{:}_highRes.dat'.format( ID ) )

indAmax = timeA.shape[0]
indBmax = timeB.shape[0]

indA = np.linspace( 0, indAmax-1, indAmax, dtype = np.int64 )
indB = np.linspace( 0, indBmax-1, indBmax, dtype = np.int64 )

fig, ax = plt.subplots( 1, 1 )
ax.set_title( r'\texttt{{{:}}}'.format( ID ) )

t = np.linspace( 5.0, min( timeA.max(), timeB.max() ), 5 )
alpha = np.linspace( 0.2, 1.0, t.shape[0] )

t = np.array( [ min( timeA.max(), timeB.max() ) ], np.float64 )
t = np.array( [ 1.0e2 ], np.float64 )
alpha = np.array( [ 1.0 ], np.float64 )

for i in range( t.shape[0] ):

    indAA = np.where( timeA < t[i] )[0]
    indBB = np.where( timeB < t[i] )[0]

    timeAA = np.copy( timeA[indAA] )
    timeBB = np.copy( timeB[indBB] )

    dataAA = np.copy( dataA[indAA] )
    dataBB = np.copy( dataB[indBB] )

    # Compute and plot FFT

    NA = timeAA.shape[0]
    NB = timeBB.shape[0]

    TA = timeAA[-1] / NA
    TB = timeBB[-1] / NB

    if( NA % 2 == 0 ):
        yA = np.abs( sp.fft.fft    ( dataAA )[1:NA//2] )**2
        xA =         sp.fft.fftfreq( NA, TA )[1:NA//2]
    else:
        yA = np.abs( sp.fft.fft    ( dataAA )[1:(NA-1)//2] )**2
        xA =         sp.fft.fftfreq( NA, TA )[1:(NA-1)//2]
    if( NA % 2 == 0 ):
        yB = np.abs( sp.fft.fft    ( dataBB )[1:NB//2] )**2
        xB =         sp.fft.fftfreq( NB, TB )[1:NB//2]
    else:
        yB = np.abs( sp.fft.fft    ( dataBB )[1:(NB-1)//2] )**2
        xB =         sp.fft.fftfreq( NB, TB )[1:(NB-1)//2]

    xA = 1.0 / xA
    xB = 1.0 / xB

    xA = xA[::-1]
    yA = yA[::-1]
    xB = xB[::-1]
    yB = yB[::-1]

    ax.plot( xA, np.abs( yA / yA.max() ), \
             'k-', \
             alpha = alpha[i], \
             label = r'(low-res) $T_{{\mathrm{{SASI}}}}$: {:.2f} ms' \
                     .format( xA[yA.argmax()] ) )

    ax.plot( xB, np.abs( yB / yB.max() ), \
             'r-', \
             alpha = alpha[i], \
             label = r'(high-res) $T_{{\mathrm{{SASI}}}}$: {:.2f} ms' \
                     .format( xB[yB.argmax()] ) )

ax.legend()
ax.set_xlabel( 'Time [ms]' )
ax.set_ylabel( 'Normalized FFT Amplitude' )

ax.grid()

plt.savefig( 'fig.{:}_FFT_Normed.png'.format( ID ), dpi = 300 )
#plt.show()
plt.close()
