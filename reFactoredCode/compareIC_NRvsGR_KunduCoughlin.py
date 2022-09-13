#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )
from scipy.optimize import curve_fit

from UtilitiesModule import GetFileArray, GetData

#### ========== User Input ==========

rootDirectory_NR \
  = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Applications/\
StandingAccretionShock_NonRelativistic/'
rootDirectory_GR \
  = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Applications/\
StandingAccretionShock_Relativistic/'

Rsh = np.array( [ '180', '150', '120', '090', '075', '060' ], str )
#Rsh = np.array( [ '020' ], str )

Mpns  = 2.8
Mdot  = '0.3'
R_PNS = 3.0
nX    = '1970'

IDs = np.array( [ '1D_M{:.1f}_Mdot{:}_Rs{:}_nX{:}' \
                  .format( Mpns, Mdot, r, nX ) for r in Rsh ] )
Rsh = np.asarray( Rsh, dtype = np.float64 )

nRows = 3
nCols = 1

yLabels = np.array( [ r'$v/v_{a}$', \
                      r'$\rho/\rho_{a}$', \
                      r'$p/\left(\rho_{a}\,v_{a}^{2}\right)$' ], str )

saveFigAs = 'fig.CompareNRvsGR_IC_KunduCoughlin.png'

verbose = False

#### ====== End of User Input =======

### Plotting

fig, axs = plt.subplots( nRows, nCols, figsize = (12,9) )
fig.suptitle( 'M{:.1f}, Mdot{:}, Rpns{:d}' \
              .format( Mpns, Mdot, np.int64( R_PNS ) ), \
              fontsize = 15, y = 0.92 )

plotFileDirectory_NR = rootDirectory_NR
plotFileDirectory_GR = rootDirectory_GR

color = [ 'red',  \
          'orange',  \
          'yellow',  \
          'green',  \
          'blue',  \
          'magenta' ]

Mpns *= 2.0e30
G_MKS = 6.673e-11
c_m  = 2.99792458e8
c_cm = 2.99792458e10

def convertLength( X1 ):
    #X = X1 * ( 1.0 + G_MKS * Mpns / ( 2.0 * c_m**2 * ( X1 * 1.0e3 ) ) )**2
    X = X1
    return X

def computeProperDistanceTo( r, rs ):
    dp = rs * ( np.sqrt( ( r * rs + r**2 ) / rs**2 ) \
                  + np.arcsinh( np.sqrt( r / rs ) ) )
    return dp

rs = G_MKS * Mpns / ( 2.0 * c_m**2  ) * 1.0e-3 # km

#Rsh   = convertLength( Rsh   )
#R_PNS = convertLength( R_PNS )

Rsh   = computeProperDistanceTo( Rsh  , rs )
R_PNS = computeProperDistanceTo( R_PNS, rs )

def fittingFunctionG( logr, m1, logr0, logrho0, m2, m3 ):
    return logrho0 + m1 * ( logr - logr0 ) \
                   + m2 * ( logr - logr0 )**2 \
                   + m3 * ( logr - logr0 )**3
def fittingFunctionN( logr, m, logr0, logrho0 ):
    return logrho0 + m * ( logr - logr0 )

for i in range( IDs.shape[0] ):

    plotFileBaseName_GR = 'GR' + IDs[i] + '.plt'
    plotFileBaseName_NR = 'NR' + IDs[i] + '.plt'

    plotFileArray_GR \
      = GetFileArray( plotFileDirectory_GR, plotFileBaseName_GR )
    plotFileArray_NR \
      = GetFileArray( plotFileDirectory_NR, plotFileBaseName_NR )

    plotFile_GR = plotFileDirectory_GR + plotFileArray_GR[-1]
    plotFile_NR = plotFileDirectory_NR + plotFileArray_NR[-1]

    time, rho_GR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_GR, 'PF_D' , verbose = verbose )
    time, rho_NR, dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_NR, 'PF_D' , verbose = verbose )
    time, p_GR  , dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_GR, 'AF_P' , verbose = verbose )
    time, p_NR  , dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_NR, 'AF_P' , verbose = verbose )
    time, v_GR  , dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_GR, 'PF_V1', verbose = verbose )
    time, v_NR  , dataUnits, X1, X2, X3, dX1, dX2, dX3, nX \
      = GetData( plotFile_NR, 'PF_V1', verbose = verbose )

    lab = IDs[i][16:21]

    ## Convert from isotropic to Schwarzchild coordinates
    eta = convertLength( X1 )

    dProper = computeProperDistanceTo( X1, rs )

    eta = np.copy( dProper )

    #eta = ( eta - R_PNS ) / ( Rsh[i] - R_PNS )

    indA = np.where( eta > Rsh[i] )[0][0]

    #ind = np.where( eta <= R_PNS )[0]
    ind = np.where( eta < Rsh[i] )[0]
    #ind = np.where( eta < 1.0e100 )[0]

    eta = np.copy( eta[ind] )

    y0_GR = np.copy( v_GR[ind,0,0] / v_GR[indA,0,0] )
    y0_NR = np.copy( v_NR[ind,0,0] / v_NR[indA,0,0] )

    rhoPrime_GR = rho_GR[indA,0,0]# * ( eta / Rsh[i] )**( -3.0 / 2.0 )
    rhoPrime_NR = rho_GR[indA,0,0]# * ( eta / Rsh[i] )**( -3.0 / 2.0 )

    y1_GR = np.copy( rho_GR[ind,0,0] / rhoPrime_GR )
    y1_NR = np.copy( rho_NR[ind,0,0] / rhoPrime_NR )

    v_GR *= 1.0e5
    v_NR *= 1.0e5

    y2_GR = np.copy( p_GR[ind,0,0] / ( rhoPrime_GR * v_GR[indA,0,0]**2 ) )
    y2_NR = np.copy( p_NR[ind,0,0] / ( rhoPrime_NR * v_NR[indA,0,0]**2 ) )

    axs[0].plot( eta, y0_GR, color = color[i], ls = '-', \
                 label = r'$v_{{a}}={:.2f}$'.format(abs(v_GR[indA,0,0])/c_cm) )
    axs[0].plot( eta, y0_NR, color = color[i], ls = '--' )
    axs[0].plot( [Rsh[i],Rsh[i]], [0.0,max(y0_GR.max(),y0_NR.max())], \
                 color = color[i] )

    axs[1].plot( eta, y1_GR, color = color[i], ls = '-' )
    axs[1].plot( eta, y1_NR, color = color[i], ls = '--' )
    axs[1].plot( [Rsh[i],Rsh[i]], [0.0,min(y1_GR.min(),y1_NR.min())], \
                 color = color[i] )

    #x = np.linspace( 10.0, Rsh[i], 1000 )
    #popt, pcov \
    #  = curve_fit( fittingFunctionG, np.log10( eta ), np.log10( y1_GR ) )
    #yy = 10**( popt[2] + popt[0] * ( np.log10( x ) - popt[1] ) \
    #                   + popt[3] * ( np.log10( x ) - popt[1] )**2 \
    #                   + popt[4] * ( np.log10( x ) - popt[1] )**3 )
    #axs[1].plot( x, yy, color = color[i], ls = '-' )
    #popt, pcov \
    #  = curve_fit( fittingFunctionN, np.log10( eta ), np.log10( y1_NR ) )
    #yy = 10**( popt[0] * ( np.log10( x ) - popt[1] ) + popt[2] )
    #popt, pcov \
    #  = curve_fit( fittingFunctionG, np.log10( eta ), np.log10( y1_NR ) )
    #yy = 10**( popt[2] + popt[0] * ( np.log10( x ) - popt[1] ) \
    #                   + popt[3] * ( np.log10( x ) - popt[1] )**2 \
    #                   + popt[4] * ( np.log10( x ) - popt[1] )**3 )
    #axs[1].plot( x, yy, color = color[i], ls = '--' )

    axs[2].plot( eta, y2_GR, color = color[i], ls = '-' )
    axs[2].plot( eta, y2_NR, color = color[i], ls = '--' )
    axs[2].plot( [Rsh[i],Rsh[i]], [0.0,min(y2_GR.min(),y2_NR.min())], \
                 color = color[i] )

xlabel = r'Proper Distance $\left[\mathrm{km}\right]$'
#xlabel = r'Radial Coordinate (Proper) $\left[\mathrm{km}\right]$'
#xlabel \
#  = r'$\left(r-R_{\textsc{pns}}\right)/\left(R_{s}-R_{\textsc{pns}}\right)$'
for i in range( nRows ):

#    axs[i].axvline( 10.0 )
    axs[i].set_ylabel( yLabels[i] )
    axs[i].axvline( R_PNS, color = 'k', ls = '--' )

    axs[i].set_xlabel( xlabel )
    axs[i].set_xscale( 'log' )

    axs[i].tick_params( axis = 'x', colors = 'k', labelsize = 15 )
    axs[i].tick_params( axis = 'y', colors = 'k', labelsize = 15 )

axs[1].set_yscale( 'log' )
axs[2].set_yscale( 'log' )

axs[0].legend( loc = 1, prop = { 'size' : 13 } )

#axs[0].set_ylim( -0.01, 0.22 )
#axs[1].set_ylim( 1.0e-3, 1.0e6 )
#axs[2].set_ylim( 1.0e-5, 1.0e7 )

#ax0ticks = np.linspace( 0.0, 0.20, 5 )
#axs[0].set_yticks( ax0ticks )
#
#ax1ticks = np.logspace( -1, 5, 4 )
#axs[1].set_yticks( ax1ticks )
#
#ax2ticks = np.logspace( -4, 6, 6 )
#axs[2].set_yticks( ax2ticks )

plt.subplots_adjust( hspace = 0.3 )
plt.savefig( saveFigAs, dpi = 300 )
#plt.show()
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
