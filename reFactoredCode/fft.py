#!/usr/bin/env python3

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

#ID = '2D_M2.8_Mdot0.3_Rs6.00e1_RPNS2.00e1'
ID = '2D_M2.8_Mdot0.3_Rs180'

ID_GR = 'GR' + ID
ID_NR = 'NR' + ID
timeGR, dataGR = np.loadtxt( 'LatFlux_{:}.dat'.format( ID_GR ) )
timeNR, dataNR = np.loadtxt( 'LatFlux_{:}.dat'.format( ID_NR ) )

indmax = min( timeGR.shape[0], timeNR.shape[0] )
ind = np.linspace( 0, indmax-1, indmax, dtype = np.int64 )
#ind = np.where( timeGR < 40.0 )[0]

timeGR = np.copy( timeGR[ind] )
timeNR = np.copy( timeNR[ind] )
dataGR = np.copy( dataGR[ind] )
dataNR = np.copy( dataNR[ind] )

dataFileNameGR = '{:}_FFT.dat'.format( ID_GR )
dataFileNameNR = '{:}_FFT.dat'.format( ID_NR )

# Compute and plot integral of data vs. time

yIntGR \
  = sp.integrate.cumulative_trapezoid \
      ( np.abs( dataGR ), x = timeGR, dx = np.diff( timeGR )[0], initial = 0.0 )
yIntNR \
  = sp.integrate.cumulative_trapezoid \
      ( np.abs( dataNR ), x = timeNR, dx = np.diff( timeNR )[0], initial = 0.0 )

fig, ax = plt.subplots( 1, 1 )
fig.suptitle( r'\texttt{{{:}}}'.format( ID ) )
ax.set_title \
  ( r'$y:=y\left(t,r\right)=\int_{S^{2}}\alpha\left(r\right)$' \
      + r'$\sqrt{\gamma}\left(r,\theta\right)$' \
      + r'$F^{1}_{S_{2}}\left(t,r,\theta\right)\,d\Omega$' )
ax.plot( timeGR, yIntGR, label = 'GR' )
ax.plot( timeNR, yIntNR, label = 'NR' )
ax.set_xlabel( 'Time [ms]' )
ax.set_ylabel( r"$\int_{0}^{t} \left|y\left(t',r\right)\right|\,dt'$" )
ax.set_yscale( 'log' )

ax.legend()
ax.grid()

#plt.savefig( 'fig.{:}_LatFlux_Integral.png'.format( ID ), dpi = 300 )
plt.show()
plt.close()

# Compute and plot FFT

NGR = timeGR.shape[0]
NNR = timeNR.shape[0]

print( "Checking Parseval's theorem" )
print( 'GR' + ID )
if( NGR % 2 == 0 ):
    yGR = np.abs( sp.fft.fft    ( dataGR )[1:NGR//2] )**2
    xGR =         sp.fft.fftfreq( NGR, 1 )[1:NGR//2]
    print( np.sum( np.abs( dataGR )**2 ) )
    print( 1/NGR * np.sum( np.abs( sp.fft.fft( dataGR )[1:] )**2 ) )
else:
    yGR = np.abs( sp.fft.fft    ( dataGR )[1:(NGR-1)//2] )**2
    xGR =         sp.fft.fftfreq( NGR, 1 )[1:(NGR-1)//2]
    print( np.sum( np.abs( dataGR )**2 ) )
    print( 1/NGR * np.sum( np.abs( sp.fft.fft( dataGR )[1:(NGR-1)//2] )**2 ) )

print( 'NR' + ID )
if( NNR % 2 == 0 ):
    yNR = np.abs( sp.fft.fft    ( dataNR )[1:NNR//2] )**2
    xNR =         sp.fft.fftfreq( NNR, 1 )[1:NNR//2]
    print( np.sum( np.abs( dataNR )**2 ) )
    print( 1/NNR * np.sum( np.abs( sp.fft.fft( dataNR )[1:] )**2 ) )
else:
    yNR = np.abs( sp.fft.fft    ( dataNR )[1:(NNR-1)//2] )**2
    xNR =         sp.fft.fftfreq( NNR, 1 )[1:(NNR-1)//2]
    print( np.sum( np.abs( dataNR )**2 ) )
    print( 1/NNR * np.sum( np.abs( sp.fft.fft( dataNR )[1:(NNR-1)//2] )**2 ) )

xGR = 1.0 / xGR
xNR = 1.0 / xNR

xGR = xGR[::-1]
yGR = yGR[::-1]
xNR = xNR[::-1]
yNR = yNR[::-1]

fig, axs = plt.subplots( 2, 1 )
fig.suptitle( r'\texttt{{{:}}}'.format( ID ) )

axs[0].plot( timeGR, dataGR, label = 'GR' )
axs[0].plot( timeNR, dataNR, label = 'NR' )

axs[1].text( 0.5, 0.6, 'GR: {:.3e} ms'.format( xGR[yGR.argmax()] ), \
             transform = axs[1].transAxes )
axs[1].text( 0.5, 0.5, 'NR: {:.3e} ms'.format( xNR[yNR.argmax()] ), \
             transform = axs[1].transAxes )

axs[1].plot( xGR, np.abs( yGR / yGR.max() ), '-' )
axs[1].plot( xNR, np.abs( yNR / yNR.max() ), '-' )

axs[0].set_xlabel( 'Time [ms]' )
axs[1].set_xlabel( 'Time [ms]' )
axs[0].set_ylabel( 'Data' )
axs[1].set_ylabel( 'Normalized FFT Amplitude' )

axs[0].legend()

axs[0].grid()
axs[1].grid()

plt.subplots_adjust( hspace = 0.3 )
#plt.savefig( 'fig.{:}_FFT_Normed.png'.format( ID ), dpi = 300 )
plt.show()
plt.close()

# Compute and plot integral of FFT vs. time

xGR = 1.0 / xGR
xNR = 1.0 / xNR

xGR = xGR[::-1]
xNR = xNR[::-1]

yIntGR \
  = sp.integrate.cumulative_trapezoid \
      ( yGR, x = xGR, dx = np.diff( xGR )[0], initial = 0.0 )
yIntNR \
  = sp.integrate.cumulative_trapezoid \
      ( yNR, x = xNR, dx = np.diff( xNR )[0], initial = 0.0 )

fig, ax = plt.subplots( 1, 1 )
fig.suptitle( r'\texttt{{{:}}}'.format( ID ) )
ax.set_title( r'$\tilde{y}\left(f\right)$' \
  + r'$:=\mathrm{FFT}\left\{y\left(t;f\right)\right\}$' )
ax.plot( xGR, yIntGR, label = 'GR' )
ax.plot( xNR, yIntNR, label = 'NR' )
ax.set_xlabel( 'Frequency [1/ms]' )
ax.set_ylabel( r"$\int_{0}^{f}\left|\tilde{y}\left(f'\right)\right|^{2}\,df'$" )
ax.set_yscale( 'log' )

ax.legend()
ax.grid()

#plt.savefig( 'fig.{:}_FFT_Integral.png'.format( ID ), dpi = 300 )
plt.show()
plt.close()
