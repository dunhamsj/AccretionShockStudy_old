#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from Relaxation import Relaxation

nX  = 640
SSi = 0
SSf = 1999
nSS = 2000

UseLogScale = True

Data = Relaxation( nX, SSi, SSf, nSS, ForceChoice = True, Overwrite = True )

ID = np.array( [ 'NR1D_M1.4_Mdot0.3_Rs180_PA0.00e-00_nX640', \
                 'GR1D_M1.4_Mdot0.3_Rs180_PA0.00e-00_nX640' ], np.str )

D = 'PF_D'
V = 'PF_V1'
P = 'AF_P'

fig, axs = plt.subplots( 3, 1 )

for i in range( ID.shape[0] ):

    IDi = ID[i]

    Time_D, Data_D = Data.GetData( D, IDi )
    Time_V, Data_V = Data.GetData( V, IDi )
    Time_P, Data_P = Data.GetData( P, IDi )

    Data.PlotRelaxationVsTime( axs[0], Time_D, Data_D, D, IDi, UseLogScale, \
                              label = IDi )
    Data.PlotRelaxationVsTime( axs[1], Time_V, Data_V, V, IDi, UseLogScale )
    Data.PlotRelaxationVsTime( axs[2], Time_P, Data_P, P, IDi, UseLogScale )

    del Time_D, Data_D
    del Time_V, Data_V
    del Time_P, Data_P

lab = 'PA0.00e-00'
fig.suptitle( 'Gaussian perturbation below shock\n{:}'.format( lab ) )

axs[0].legend(loc=1)
axs[0].set_ylim( 1.0e-6, 5.0 )
axs[1].set_ylim( 1.0e-6, 5.0 )
axs[2].set_ylim( 1.0e-6, 5.0 )

#plt.show()

SaveFigAs = 'fig.Relaxation_{:}.png'.format( lab )
plt.savefig( SaveFigAs, dpi = 300 )

plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
