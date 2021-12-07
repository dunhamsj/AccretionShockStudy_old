#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

DataShape = (602,320,64)

Data = np.loadtxt( '.GR2D_M1.4_Mdot0.3_Rs180_PA0.000_nX320x064_DF_TCI.dat' ).reshape( DataShape )

D = Data[:,:,32]

for i in range( D.shape[0] ):
    for j in range( D.shape[1] ):
        D[i,j] = max( D[i,j], 1.0e-2 )
D = np.log10( D )

fig, ax = plt.subplots()

extent = [ 40.0, 360.0, 0.0, 600.0 ]

im = ax.imshow( D, \
                extent = extent )

ax.set_xlabel( 'Radial Coordinate [km]' )
ax.set_ylabel( 'Time [ms]' )
cbar = fig.colorbar( im )
cbar.set_label( 'Log10( TCI )' )

#plt.show()

plt.savefig( 'fig.HeatMap_TCI.png', dpi = 300 )
