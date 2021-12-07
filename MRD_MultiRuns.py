#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from MaximumRelativeDifference import MaxRelDiff

nX  = 640
SSi = 0
SSf = 1999
nSS = 2000

UseLogScale = True

MRD = MaxRelDiff( nX, SSi, SSf, nSS, ForceOverwrite = True, Overwrite = False )

ID_0 = 'NR1D_M1.4_Mdot0.3_Rs180_PA0.00e-00_nX640'
ID_1 = 'NR1D_M1.4_Mdot0.3_Rs180_PA1.00e-01_nX640'
ID_2 = 'NR1D_M1.4_Mdot0.3_Rs180_PA1.00e-03_nX640'
ID_3 = 'NR1D_M1.4_Mdot0.3_Rs180_PA1.00e-05_nX640'

D = 'PF_D'
Time_D_0, Data_D_0 = MRD.GetData( D, ID_0 )
Time_D_1, Data_D_1 = MRD.GetData( D, ID_1 )
Time_D_2, Data_D_2 = MRD.GetData( D, ID_2 )
Time_D_3, Data_D_3 = MRD.GetData( D, ID_3 )

V = 'PF_V1'
Time_V_0, Data_V_0 = MRD.GetData( V, ID_0 )
Time_V_1, Data_V_1 = MRD.GetData( V, ID_1 )
Time_V_2, Data_V_2 = MRD.GetData( V, ID_2 )
Time_V_3, Data_V_3 = MRD.GetData( V, ID_3 )

P = 'AF_P'
Time_P_0, Data_P_0 = MRD.GetData( P, ID_0 )
Time_P_1, Data_P_1 = MRD.GetData( P, ID_1 )
Time_P_2, Data_P_2 = MRD.GetData( P, ID_2 )
Time_P_3, Data_P_3 = MRD.GetData( P, ID_3 )

fig, axs = plt.subplots( 3, 1 )

lab = 'NR1D_M1.4_Mdot0.3_Rs180_nX640'
fig.suptitle( 'Gaussian perturbation below shock\n{:}'.format( lab ) )

MRD.PlotMaxRelDiffVsTime( axs[0], Time_D_0, Data_D_0, D, ID_0, \
                          UseLogScale, label = 'PA0.00e-00' )
MRD.PlotMaxRelDiffVsTime( axs[0], Time_D_1, Data_D_1, D, ID_1, \
                          UseLogScale, label = 'PA1.00e-01' )
MRD.PlotMaxRelDiffVsTime( axs[0], Time_D_1, Data_D_2, D, ID_2, \
                          UseLogScale, label = 'PA1.00e-03' )
MRD.PlotMaxRelDiffVsTime( axs[0], Time_D_1, Data_D_3, D, ID_3, \
                          UseLogScale, label = 'PA1.00e-05' )

MRD.PlotMaxRelDiffVsTime( axs[1], Time_V_0, Data_V_0, V, ID_0, UseLogScale )
MRD.PlotMaxRelDiffVsTime( axs[1], Time_V_1, Data_V_1, V, ID_1, UseLogScale )
MRD.PlotMaxRelDiffVsTime( axs[1], Time_V_2, Data_V_2, V, ID_2, UseLogScale )
MRD.PlotMaxRelDiffVsTime( axs[1], Time_V_3, Data_V_3, V, ID_3, UseLogScale )

MRD.PlotMaxRelDiffVsTime( axs[2], Time_P_0, Data_P_0, P, ID_0, UseLogScale )
MRD.PlotMaxRelDiffVsTime( axs[2], Time_P_1, Data_P_1, P, ID_1, UseLogScale )
MRD.PlotMaxRelDiffVsTime( axs[2], Time_P_2, Data_P_2, P, ID_2, UseLogScale )
MRD.PlotMaxRelDiffVsTime( axs[2], Time_P_3, Data_P_3, P, ID_3, UseLogScale )

axs[0].legend()
axs[0].set_ylim( 1.0e-6, 5.0 )
axs[1].set_ylim( 1.0e-6, 5.0 )
axs[2].set_ylim( 1.0e-6, 5.0 )

#plt.show()

SaveFigAs = 'fig.MaxRelDiff_{:}.png'.format( lab )
plt.savefig( SaveFigAs, dpi = 300 )

plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
