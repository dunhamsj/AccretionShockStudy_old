#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

Field = 'DivV2'

#a = np.arange( 0.2, 1.0+0.05, 0.05 )
#a = np.arange( 0.65, 0.90, 0.05 )
a = np.arange( 0.80, 0.91, 0.01 )

N = a.shape[0]

annuli = []
for n in range( N - 1 ):

#    annuli.append( '{:.2f}-{:.2f}'.format( a[n], a[n+1] ) )
    annuli.append( '{:.2f}-{:.2f}'.format( a[n], 0.9 ) )

Root = '.GR2D_M1.4_Mdot0.3_Rs180_PA0.004_nX320x064_{:s}_'.format( Field )
suffix = '_PowersInLegendreModes.dat'

nRows = np.int64( np.sqrt( np.float64( len( annuli ) ) ) )

if np.sqrt( np.float64( len( annuli ) ) ) > nRows: nRows += 1
nCols = nRows

fig, axs = plt.subplots( nRows, nCols, figsize = (16,9) )
tit = r'Power in $\ell=1$ mode vs. time for power computed in different annuli'
tit2 = '\nusing {:s} as variable'.format( Field )
fig.suptitle( tit + tit2 )

k = -1

for i in range( axs.shape[0] ):

    for j in range( axs.shape[1] ):

        k += 1

        if k > N - 2: continue

        t, P0, P1, P2 = np.loadtxt( Root + annuli[k] + suffix )

        ind = np.where( t > 10.0 )[0]

        axs[i,j].semilogy( t[ind], P1[ind], label = annuli[k] )

        axs[i,j].legend()
#plt.savefig( 'fig.PowerInAnnuli_{:}.png'.format( Field ), dpi = 300 )
plt.show()
