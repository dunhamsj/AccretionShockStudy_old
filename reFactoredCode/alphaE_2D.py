#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )
from scipy.optimize import curve_fit

#from MakeDataFile import MakeDataFile, ReadHeader

rootDirectory = '/lump/data/accretionShockStudy/newRuns/newProductionRuns/'

ID_1D = 'GR1D_M1.4_Rpns040_Rs180_Mdot0.3'
ID_2D = 'GR2D_M1.4_Rpns040_Rs180_Mdot0.3'

plotfileDirectory_1D = rootDirectory + ID_1D + '/'
plotfileDirectory_2D = rootDirectory + ID_2D + '/'

tallyFile_1D = plotfileDirectory_1D + ID_1D + '.Tally_Energy.dat'
tallyFile_2D = plotfileDirectory_2D + ID_2D + '.Tally_Energy.dat'

t_1D, Eint_1D, EOG_1D, Einit_1D, dE_1D \
  = np.loadtxt( tallyFile_1D, skiprows = 1, unpack = True )
t_2D, Eint_2D, EOG_2D, Einit_2D, dE_2D \
  = np.loadtxt( tallyFile_2D, skiprows = 1, unpack = True )

ind_1D = np.where( ( t_1D <= 6.0e2 ) & ( t_1D >= 0.0 ) )[0]
ind_2D = np.where( ( t_2D <= 6.0e2 ) & ( t_2D >= 0.0 ) )[0]

t_1D    = np.copy( t_1D   [ind_1D] )
Eint_1D = np.copy( Eint_1D[ind_1D] )
dE_1D   = np.copy( dE_1D  [ind_1D] )

t_2D    = np.copy( t_2D   [ind_2D] )
Eint_2D = np.copy( Eint_2D[ind_2D] )
dE_2D   = np.copy( dE_2D  [ind_2D] )

dr   = 0.5
Eend_1D = dE_1D[-1] / Eint_1D[0]
Eend_2D = dE_2D[-1] / Eint_2D[0]

fig, ax = plt.subplots( 1, 1 )
ax.set_title( r'$\texttt{{{:}}}$'.format( ID_2D[5:] ) )
ax.plot( t_1D, dE_1D / Eint_1D[0], '-', label = '1D' )
ax.plot( t_2D, dE_2D / Eint_2D[0], '-', label = '2D' )
ax.set_xlabel( r'$t\ \left[\mathrm{ms}\right]$' )
ax.set_ylabel( \
r'$\left(y\left(t\right)-y\left(0\right)+\widehat{F}_{y_{\textrm{off-grid}}}\right)/y\left(0\right)$', fontsize = 12 )

ax.text( 0.1, 0.7, r'$y\left(t\right):=\int_{V}\alpha\,E\left(t\right)\,dV$', \
         transform = ax.transAxes, fontsize = 15 )
ax.legend()
#plt.show()
plt.savefig( '/home/kkadoogan/fig.IntegralOfAlphaEvsTime_2D.png', dpi = 300 )

import os
os.system( 'rm -rf __pycache__ ' )
