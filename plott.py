#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt( 'TMPdata.dat' )
for i in range( 1, 5 ):
    plt.semilogy( data[0], data[i] )
plt.show()
