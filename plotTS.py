#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'Publication.sty' )

TS = np.loadtxt( 'TS.dat', delimiter = ',' )

Gm = TS[:,0]
tauAd = TS[:,1]
tauAc = TS[:,2]
tauSASI = tauAd + tauAc

plt.title( 'NR1D_M2.0_Mdot0.3_Rs150' )
#plt.semilogy( Gm, tauAd, '.', label = r'$\tau_{\mathrm{Advective}}$' )
plt.semilogy( Gm, tauAd/tauAc, '.', label = r'$\tau_{\mathrm{Acoustic}}$' )
#plt.semilogy( Gm, tauSASI, '.', label = r'$\tau_{\mathrm{SASI}}$' )
plt.grid()
plt.xlabel( 'Gamma' )
plt.ylabel( 'Time [ms]' )
plt.legend()
#plt.savefig( 'fig.VaryGm.png', dpi = 300 )
plt.show()
