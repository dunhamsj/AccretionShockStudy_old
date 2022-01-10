#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from MaximumRelativeDifference import MaxRelDiff

nX  = 640
SSi = 0
SSf = 1999
nSS = 2000

UseLogScale = True

MRD = MaxRelDiff( nX, SSi, SSf, nSS, ForceChoice = True, Overwrite = True )

PA = np.array( [ 'PA0.00e-00', \
                 'PA1.00e-01', \
                 'PA1.00e-03', \
                 'PA1.00e-04', \
                 'PA1.00e-05' ], np.str )

D = 'PF_D'
V = 'PF_V1'
P = 'AF_P'

fig, axs = plt.subplots( 3, 1 )

for i in range( PA.shape[0] ):

    ID = 'NR1D_M1.4_Mdot0.3_Rs180_{:}_nX640'.format( PA[i] )

    Time_D, Data_D = MRD.GetData( D, ID )
    Time_V, Data_V = MRD.GetData( V, ID )
    Time_P, Data_P = MRD.GetData( P, ID )

    MRD.PlotMaxRelDiffVsTime( axs[0], Time_D, Data_D, D, ID, UseLogScale, \
                              label = PA[i] )
    MRD.PlotMaxRelDiffVsTime( axs[1], Time_V, Data_V, V, ID, UseLogScale )
    MRD.PlotMaxRelDiffVsTime( axs[2], Time_P, Data_P, P, ID, UseLogScale )

    del Time_D, Data_D
    del Time_V, Data_V
    del Time_P, Data_P

lab = 'NR1D_M1.4_Mdot0.3_Rs180_nX640'
fig.suptitle( 'Gaussian perturbation below shock\n{:}'.format( lab ) )

axs[0].legend(loc=1)
axs[0].set_ylim( 1.0e-6, 5.0 )
axs[1].set_ylim( 1.0e-6, 5.0 )
axs[2].set_ylim( 1.0e-6, 5.0 )

#plt.show()

SaveFigAs = 'fig.MaxRelDiff_{:}.png'.format( lab )
plt.savefig( SaveFigAs, dpi = 300 )

plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
