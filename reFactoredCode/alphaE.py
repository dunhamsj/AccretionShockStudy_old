#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )
from scipy.optimize import curve_fit

#from MakeDataFile import MakeDataFile, ReadHeader

rootDirectory = '/lump/data/accretionShockStudy/newRuns/alphaE_resolutionStudy/'

ID = 'GR1D_M1.4_Rpns040_Rs150_Mdot0.3'

suffix = [ '.nX0320', '.nX0640', '.nX1280' ]

#alphaE = 'alphaE'
#SqrtGm = 'GF_SqrtGm'
#
#SSi = -1
#SSf = -1
#nSS = -1
#
#fcD = True
#fcF = True

N = len( suffix )

t    = np.empty( N, object )
Eint = np.empty( N, object )
dE   = np.empty( N, object )

dr   = np.empty( N, np.float64 )
Eend = np.empty( N, np.float64 )

for i in range( N ):

    IDD = ID + suffix[i]

    plotfileDirectory = rootDirectory + IDD + '/'

#    dataDirectory = '.{:}/'.format( IDD )
#
#    plotfileBaseName = IDD + '.plt'
#
#    plotfileArray \
#      = MakeDataFile( alphaE, plotfileDirectory, dataDirectory, \
#                      plotfileBaseName, 'spherical', \
#                      SSi = SSi, SSf = SSf, nSS = nSS, \
#                      UsePhysicalUnits = True, \
#                      MaxLevel = -1, \
#                      forceChoiceD = fcD, owD = False, \
#                      forceChoiceF = fcF, owF = False, \
#                      Verbose = True )
#    plotfileArray \
#      = MakeDataFile( SqrtGm, plotfileDirectory, dataDirectory, \
#                      plotfileBaseName, 'spherical', \
#                      SSi = SSi, SSf = SSf, nSS = nSS, \
#                      UsePhysicalUnits = True, \
#                      MaxLevel = -1, \
#                      forceChoiceD = True, owD = False, \
#                      forceChoiceF = fcF , owF = False, \
#                      Verbose = True )
#
#    plotfileArray = np.copy( plotfileArray[:-1] )
#    nSS = plotfileArray.shape[0]
#
#    time     = np.zeros( nSS, np.float64 )
#    integral = np.zeros( nSS, np.float64 )
#
#    dX2 = np.pi
#    dX3 = 2.0 * np.pi
#    for iSS in range( nSS ):
#
#        fileDirectory = dataDirectory + plotfileArray[iSS] + '/'
#
#        timeFile   = fileDirectory + '{:}.dat'.format( 'Time' )
#        X1File     = fileDirectory + '{:}.dat'.format( 'X1' )
#        dX1File    = fileDirectory + '{:}.dat'.format( 'dX1' )
#
#        time[iSS] = np.loadtxt( timeFile )
#        X1_C      = np.loadtxt( X1File   )
#        dX1       = np.loadtxt( dX1File  )
#
#        alphaEFile = fileDirectory + '{:}.dat'.format( alphaE )
#        SqrtGmFile = fileDirectory + '{:}.dat'.format( SqrtGm )
#
#        alphaEdata = np.loadtxt( alphaEFile )
#        SqrtGmdata = np.loadtxt( SqrtGmFile )
#
#        alphaEShape, alphaEUnits, alphaEminVal, alphaEmaxVal \
#          = ReadHeader( alphaEFile )
#        SqrtGmShape, SqrtGmUnits, SqrtGmminVal, SqrtGmmaxVal \
#          = ReadHeader( SqrtGmFile )
#
#        integral[iSS] = np.sum( alphaEdata * SqrtGmdata * dX1 ) * dX2 * dX3

    tallyFile = plotfileDirectory + IDD + '.Tally_Energy.dat'
    tt, Eintt, EOG, Einit, dEt \
      = np.loadtxt( tallyFile, skiprows = 1, unpack = True )

    ind = np.where( ( tt <= 1.0e3 ) & ( tt >= 0.0 ) )[0]

    print(ind.shape )
    t   [i] = np.copy( tt   [ind] )
    Eint[i] = np.copy( Eintt[ind] )
    dE  [i] = np.copy( dEt  [ind] )

    dr  [i] = ( 3.60e2 - 4.00e1 ) / np.int64( suffix[i][3:] )
    Eend[i] = dE[i][-1] / Eint[i][0]

fig, axs = plt.subplots( 2, 1 )
fig .suptitle( r'$\texttt{{{:}}}$'.format( ID ) )
for i in range( N ):
    axs[0].plot( t[i], Eint[i] / Eint[i][0], '-', label = suffix[i][1:] )
    axs[1].plot( t[i], dE  [i] / Eint[i][0], '-' )
axs[0].set_title( r'$y:=\int_{V}\alpha\,E\,dV$' )
axs[0].set_ylabel( r'$y\left(t\right)/y\left(0\right)$' )
axs[0].set_xticklabels('')
axs[0].legend()
axs[1].set_xlabel( r'$t\ \left[\mathrm{ms}\right]$' )
axs[1].set_ylabel( \
r'$\left(y\left(t\right)-y\left(0\right)+\widehat{F}_{y_{\textrm{off-grid}}}\right)/y\left(0\right)$', fontsize = 12 )
fig.subplots_adjust( hspace = 0.0 )

#plt.show()
plt.savefig( '/home/kkadoogan/fig.IntegralOfAlphaEvsTime.png', dpi = 300 )

plt.close()

fig, ax = plt.subplots( 1, 1 )
fig .suptitle( r'$\texttt{{{:}}}$'.format( ID ) )
Eend = Eend / Eend.max()
ax.plot( dr, Eend, 'o' )
ax.set_title( r'$y:=\int_{V}\alpha\,E\,dV$' )
ax.set_ylabel( \
r'$\left(y\left(t=1000\,\mathrm{ms}\right)-y\left(t=0\right)+\widehat{F}_{y_{\textrm{off-grid}}}\right)/y\left(t=0\right)$', fontsize = 12 )
ax.set_xlabel( r'$dr\ \left[\mathrm{km}\right]$' )
def y( x, m, b ):
    return m * x + b
popt, pcov = curve_fit( y, dr, Eend )

m = popt[0]
b = popt[1]

ax.plot( dr, m * dr + b, label = r'$m={:.2e}$'.format( m ) )

#ax.plot( dr, dr )
#ax.set_xscale( 'log' )
#ax.set_yscale( 'log' )
ax.legend()

plt.show()
#plt.savefig( '/home/kkadoogan/fig.ConvergenceRate.png', dpi = 300 )

plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
