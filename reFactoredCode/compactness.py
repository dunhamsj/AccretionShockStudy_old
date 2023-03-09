#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )
from scipy.optimize import curve_fit
from os.path import isfile

# Compute compactness

usePeriods = True

R0 = 10.0 # km

M  = np.array( [ 1.4, 2.0, 2.8  ], np.float64 )
R  = np.array( [ 20.0, 40.0 ], np.float64 )
Rs = np.array( [ 60, 75, 90, 120, 150, 180 ], np.int64 )

xi = np.empty( (M.shape[0],R.shape[0]), np.float64 )
for m in range( M.shape[0] ):
    for r in range( R.shape[0] ):
        xi[m,r] = M[m] / ( R[r] / R0 )

# Read in growth rates

growthRateFileDirectory \
  = '/home/kkadoogan/Work/accretionShockPaper/plottingData/'

def remove_sign( value ):
    s = '{:.2e}'.format( value ).replace( 'e+0', 'e' )
    return s

def fit( x, m, b ):
    return m * x + b

y = np.zeros( (M.shape[0],Rs.shape[0]), np.float64 )

for r in range( R.shape[0] ):

    for rs in range( Rs.shape[0]-1, -1, -1 ):

        if( ( r == 0 ) & ( rs > 2 ) ): continue
        if( ( r == 1 ) & ( rs < 3 ) ): continue

        for m in range( M.shape[0] ):

            if( ( r == 0 ) & ( m < 2 ) ): continue

            if rs < 3:
                model = 'M{:.1f}_Mdot0.3_Rs{:}_RPNS2.00e1' \
                        .format( M[m], remove_sign( Rs[rs] ) )
            else:
                model = 'M{:.1f}_Mdot0.3_Rs{:d}'.format( M[m], Rs[rs] )

            file = '{:}NR2D_{:}_Fit.dat' \
                   .format( growthRateFileDirectory, model )

            data_NR \
              = np.loadtxt( '{:}NR2D_{:}_Fit.dat' \
                            .format( growthRateFileDirectory, model ) )
            data_GR \
              = np.loadtxt( '{:}GR2D_{:}_Fit.dat' \
                            .format( growthRateFileDirectory, model ) )

            if rs < 3:
                model = model[0:4] + model[12:19]
            else:
                model = model[0:4] + model[-6:]

            if usePeriods:

                G_NR = 2.0 * np.pi / data_NR[4]
                G_GR = 2.0 * np.pi / data_GR[4]

            else:

              G_NR = data_NR[3]
              G_GR = data_GR[3]

            y[m,rs] = G_GR - G_NR

fig, ax = plt.subplots( 1, 1 )

# colorblind-friendly palette: https://gist.github.com/thriveth/8560036
color = ['#377eb8', '#ff7f00', '#4daf4a', \
         '#f781bf', '#a65628', '#984ea3', \
         '#999999', '#e41a1c', '#dede00']
alpha = np.linspace( 0.2, 1.0, M.shape[0] )

s = [ '*', '.' ]
for r in range( R.shape[0] ):

    if r == 0: continue

    for rs in range( Rs.shape[0]-1, -1, -1 ):

        if( ( r == 0 ) & ( rs > 2 ) ): continue
        if( ( r == 1 ) & ( rs < 3 ) ): continue

        for m in range( M.shape[0] ):

            if( ( r == 0 ) & ( m < 2 ) ): continue

            if rs < 3:
                model = 'M{:.1f}_Mdot0.3_Rs{:}_RPNS2.00e1' \
                        .format( M[m], remove_sign( Rs[rs] ) )
            else:
                model = 'M{:.1f}_Mdot0.3_Rs{:d}'.format( M[m], Rs[rs] )

            file = '{:}NR2D_{:}_Fit.dat' \
                   .format( growthRateFileDirectory, model )

            if rs < 3:
                model = model[0:4] + '_Rs{:}_Rpns020' \
                        .format( str( Rs[rs] ).zfill( 3 ) )
            else:
                model = model[0:4] + model[-6:] + '_Rpns040'

            ax.plot( xi[m,r], y[m,rs], s[r], color = color[rs%3], \
                     alpha = alpha[m], \
                     label = r'\texttt{{{:}}}'.format( model ) )

#        popt, pcov = curve_fit( fit, xi[:,r], y[:,rs] )
#        ax.plot( xi[:,r], popt[0] * xi[:,r] + popt[1], color = color[rs%3], \
#                 label = r'$m={:.2f}$'.format( popt[0] ) )

if usePeriods:
    loc = (0.1,0.3)
    ax.set_ylabel( r'$T_{\mathrm{GR}}-T_{\mathrm{NR}}\ \left[\mathrm{ms}\right]$' )
else:
    loc = (0.1,0.1)
    ax.set_ylabel( r'$\omega_{\mathrm{GR}}-\omega_{\mathrm{NR}}\ \left[\mathrm{ms}^{-1}\right]$' )

ax.legend( loc = loc )

ax.set_xlabel( r'$\xi:=M/M_{\odot}/\left(R_{\mathrm{PNS}}/10\,\mathrm{km}\right)$' )

plt.show()

#if usePeriods:
#    plt.savefig( '/home/kkadoogan/fig.Compactness_AllRuns_Period.png', dpi = 300 )
#else:
#    plt.savefig( '/home/kkadoogan/fig.Compactness_AllRuns_GrowthRate.png', dpi = 300 )
