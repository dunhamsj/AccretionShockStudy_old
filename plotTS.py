#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'Publication.sty' )

TS = np.loadtxt( 'TS.dat', delimiter = ',' )

Gamma = TS[::-1,0]
tauAd = TS[::-1,1]
tauAc = TS[::-1,2]
tauSASI = tauAd + tauAc

fig, ax = plt.subplots()

ax.set_title( 'NR1D_M2.0_Mdot0.3_Rs150' )
#ax.semilogx( tauAd, Gamma, '.', label = r'$\tau_{\mathrm{Advective}}$' )
#ax.semilogx( tauAc, Gamma, '.', label = r'$\tau_{\mathrm{Acoustic}}$' )

Ratio = tauAd / tauAc

ax.semilogx( Ratio, Gamma, '.-', lw = 2, label = r'$\tau_{Ad}/\tau_{Ac}$' )

def Fit( desiredValue ):
  ind = np.where( Ratio < desiredValue )[0][-1]
  logRatio1 = np.log10( Ratio[ind  ] )
  logRatio2 = np.log10( Ratio[ind+1] )
  Gamma1 = Gamma[ind  ]
  Gamma2 = Gamma[ind+1]
  m = ( Gamma2 - Gamma1 ) / ( logRatio2 - logRatio1 )
  return Gamma1, m, logRatio1

def interpolate( logRatio, Gamma1, m, logRatio1 ):
  return Gamma1 + m * ( logRatio - logRatio1 )

# ratio = 10
ax.axvline( 1.0e1, c = 'k' )
Gamma1, m, logRatio1 = Fit( 1.0e1 )
logRatio = np.linspace( 0.9, 1.1, 1000 )
ax.plot( 10**( logRatio ), interpolate( logRatio, Gamma1, m, logRatio1 ), \
         'm' )
Gam = interpolate( 1.0, Gamma1, m, logRatio1 )
ax.text( 1.0e1, Gam, '{:.3f}'.format( Gam ) )
print( 'Ratio = 1.0e1, Gamma: {:.16e}'.format( Gam ) )

# ratio = 100
ax.axvline( 1.0e2, c = 'k' )
Gamma1, m, logRatio1 = Fit( 1.0e2 )
logRatio = np.linspace( 1.9, 2.1, 1000 )
ax.plot( 10**( logRatio ), interpolate( logRatio, Gamma1, m, logRatio1 ), \
         'm' )
Gam = interpolate( 2.0, Gamma1, m, logRatio1 )
ax.text( 1.0e2, Gam, '{:.3f}'.format( Gam ) )
print( 'Ratio = 1.0e2, Gamma: {:.16e}'.format( Gam ) )

# ratio = 1000
ax.axvline( 1.0e3, c = 'k' )
Gamma1, m, logRatio1 = Fit( 1.0e3 )
logRatio = np.linspace( 2.9, 3.1, 1000 )
ax.plot( 10**( logRatio ), interpolate( logRatio, Gamma1, m, logRatio1 ), \
         'm' )
Gam = interpolate( 3.0, Gamma1, m, logRatio1 )
ax.text( 1.0e3, Gam, '{:.3f}'.format( Gam ) )
print( 'Ratio = 1.0e3, Gamma: {:.16e}'.format( Gam ) )

plt.grid()
plt.xlim( 1.0e0, 2.0e3 )
plt.ylabel( 'Gamma' )
plt.xlabel( 'Time [ms]' )
plt.legend()
plt.savefig( 'fig.VaryGamma.png', dpi = 300 )
#plt.show()
