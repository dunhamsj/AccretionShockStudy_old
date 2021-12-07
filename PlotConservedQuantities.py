#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

DataDirectory = '/Users/dunhamsj/Research/Data/AccretionShockParameterStudy/'
ID = 'GR2D_M1.4_Mdot0.3_Rs180_PA0.040_nX320x064_smooth'

DataDirectory += ID + '/'

DataMass   = np.loadtxt( DataDirectory + ID + '.plt_StandingAccretionShock_Relativis_Tally_BaryonicMass.dat', skiprows = 1 )
DataEnergy = np.loadtxt( DataDirectory + ID + '.plt_StandingAccretionShock_Relativis_Tally_Energy.dat', skiprows = 1 )

Time = DataMass[:,0]

Mass   = ( DataMass  [:,1] + DataMass  [:,2] ) / DataMass  [:,3]
Energy = ( DataEnergy[:,1] + DataEnergy[:,2] ) / DataEnergy[:,3]

fig, ax = plt.subplots( 1, 1 )
fig.suptitle( 'Relative to initial' )

ax.plot( Time, Mass, label = 'Baryonic Mass' )
ax.plot( Time, Energy, label = 'Energy' )
ax.set_xlabel( 'Time [ms]' )
ax.legend()

plt.show()
