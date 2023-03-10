#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt( '/lump/data/accretionShockStudy/newRuns/newProductionRuns/NR2D_M1.4_Rpns040_Rs180_Mdot0.3/NR1D_M1.4_Rpns040_Rs180_Mdot0.3.IC', skiprows = 3 )

D = data[0::3].flatten()
V = data[1::3].flatten()
P = data[2::3].flatten()
plt.plot( V )
plt.show()
