#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from FitDataToModel import FittingFunction

ind0 = 1
ind1 = 200

fit \
  = np.loadtxt \
    ( '.NR2D_M2.8_Mdot0.3_Rs6.00e1_RPNS2.00e1_Fit.dat' )
data \
  = np.loadtxt \
    ( '.NR2D_M2.8_Mdot0.3_Rs6.00e1_RPNS2.00e1_LegendrePowerSpectrum.dat' )
data2 \
  = np.loadtxt \
    ( '.NR2D_M2.8_Mdot0.3_Rs6.00e1_RPNS2.00e1_LegendrePowerSpectrum.dat2.dat' )

t0     = fit[0]
t1     = fit[1]
logF1  = fit[2]
omegaR = fit[3]
omegaI = fit[4]
delta  = fit[5]

tD = data[0][ind0:ind1]
tF = np.linspace( t0, t1, tD.shape[0] )

P0 = data[1][ind0:ind1]
P1 = data[2][ind0:ind1]
P2 = data[3][ind0:ind1]
P3 = data[4][ind0:ind1]
P4 = data[5][ind0:ind1]
yF = FittingFunction( tF, logF1, omegaR, omegaI, delta )

plt.semilogy( tD, P0, '.', label = 'P0' )
plt.semilogy( tD, P1, '.', label = 'P1' )
plt.semilogy( tD, P2, '.', label = 'P2' )
plt.semilogy( tD, P3, '.', label = 'P3' )
plt.semilogy( tD, P4, '.', label = 'P4' )
plt.semilogy( data2[0][ind0:ind1], data2[1][ind0:ind1]**2, '-', label = 'data' )
#plt.plot( tF, yF, '-' )
#plt.yscale( 'symlog', linthresh = 1.0e65 )
plt.legend()
plt.show()
