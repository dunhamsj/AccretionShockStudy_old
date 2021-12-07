#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def SetUnits( UsePhysicalUnits = True ):

    if UsePhysicalUnits:

        c = 2.99792458e8
        G = 6.673e-11

    else:

        c = 1.0
        G = 1.0

    Meter    = 1.0
    Second   = c
    Kilogram = G / c**2
    Gram = 1.0e-3 * Kilogram
    Centimeter = 1.0e-2 * Meter
    Kilometer = 1.0e3 * Meter

    return Meter, Second, Kilogram, Gram, Centimeter, Kilometer

def CreateMesh( xL, xU, nX, nN ):

    if nN == 3:

        xq \
          = np.array( [ -np.sqrt( 3.0 / 20.0 ), 0.0, +np.sqrt( 3.0 / 20.0 ) ], \
                      np.float64 )

    dX = ( xU - xL ) / np.float64( nX )

    X1C = np.linspace( xL[0] + 0.5 * dX[0], xU[0] - 0.5 * dX[0], nX[0] )
    X2C = np.linspace( xL[1] + 0.5 * dX[1], xU[1] - 0.5 * dX[1], nX[1] )
    X3C = np.linspace( xL[2] + 0.5 * dX[2], xU[2] - 0.5 * dX[2], nX[2] )

    X1 = np.empty( (nX[0]*nN), np.float64 )
    X2 = np.empty( (nX[1]*nN), np.float64 )
    X3 = np.empty( (nX[2]*nN), np.float64 )

    for i in range( nX[0] ):
        ind = np.linspace( nN*i, nN*i+nN-1, nN, dtype = np.int64 )
        X1[ind] = xq * dX[0] + X1C[i]

    for i in range( nX[1] ):
        ind = np.linspace( nN*i, nN*i+nN-1, nN, dtype = np.int64 )
        X2[ind] = xq * dX[1] + X2C[i]

    for i in range( nX[2] ):
        ind = np.linspace( nN*i, nN*i+nN-1, nN, dtype = np.int64 )
        X3[ind] = xq * dX[2] + X3C[i]

    return X1C, X2C, X3C, X1, X2, X3, dX

Meter, Second, Kilogram, Gram, Centimeter, Kilometer \
  = SetUnits( UsePhysicalUnits = True )

xL = np.array( [ 40.0, 0.0, 0.0 ], np.float64 )
xU = np.array( [ 360.0, np.pi, 2.0 * np.pi ], np.float64 )

nX = np.array( [ 320, 64, 1 ], np.int64 )
nN = 3
nDOFX = nX * nN

X1C, X2C, X3C, X1, X2, X3, dX = CreateMesh( xL, xU, nX, nN )

# Get 1D data

Data = np.loadtxt( 'GR1D_M1.4_Mdot0.3_Rs180_PA0.000_nX320.IC' )

rho1D = np.empty( (nX[0],nN), np.float64 )
v1D   = np.empty( (nX[0],nN), np.float64 )

for i in range( nX[0] ):

    rho1D[i] = Data[3+3*i] / ( Gram / Centimeter**3 )
    v1D  [i] = Data[3+3*i+1] / ( Kilometer / Second )

rho0 = rho1D.flatten()
v0   = v1D.flatten()

# Create grid of perturbation amplitudes and widths

rho1  = np.array( [ 0.01, 0.1, 0.5 ], np.float64 )
sigma = np.array( [ 1.0, 2.0, 4.0, 8.0, 16.0 ], np.float64 )

# Central location of perturbation

r0 = 260.0

Rs = 180.0

ind = np.where( X1 > Rs )[0]
indR = ind[0]

# Plotting (1D)

#fig, axs = plt.subplots( sigma.shape[0], rho1.shape[0], figsize = (16,9) )
#
#title  = r'$\rho(r)=\rho_0(r)+\rho_1\,e^{-(r-r_0)^2/(2\sigma^2)}$'
#title += '\n'
#title += r'$r_0=$ {:d} km'.format( np.int64( r0 ) )
#
#fig.suptitle( title, fontsize = 15 )
#
#xT, yT = 0.6, 0.7
#
#for i in range( sigma.shape[0] ):
#
#    for j in range( rho1.shape[0] ):
#
#        rho11 =  rho1[j] * rho0[ind[0]]
#
#        ind = np.where(   ( X1 < r0 + 0.5 * sigma[i] ) \
#                        & ( X1 > r0 - 0.5 * sigma[i] ) )[0]
#        drho = np.zeros( (nDOFX[0]), np.float64 )
#        drho[ind] = rho11
#        axs[i,j].semilogy( X1, rho0 + drho )
#
#        drho = rho11 * np.exp( -( X1 - r0 )**2 / ( 2.0 * sigma[i]**2 ) )
#        axs[i,j].semilogy( X1, rho0 + drho )
#
#        print( '\n' )
#        print( 'sigma = {:.2f}'.format( sigma[i] ) )
#        print( 'rho1  = {:.2f}'.format( rho1[j] ) )
#        print( 'drho0(Rs)/rho0(Rs) = {:.3e}'.format \
#               ( drho[indR] / rho0[indR] ) )
#
#        text = r'$\sigma$ = {:.0f} km'.format( sigma[i] )
#        text += '\n'
#        text += r'$\rho_1$ = {:.2f} $\rho_0(r_0)$'.format( rho1[j] )
#        axs[i,j].text( xT, yT, text, transform = axs[i,j].transAxes )
#
##        axs[i,j].set_ylim( 1.0e7, 3.0e7 )
##        axs[i,j].set_xlim( 1.8e2, 3.6e2 )
#        if j == 0:
#            axs[i,j].set_ylabel( r'$\rho$ [g/cc]' )
#        else:
#            axs[i,j].yaxis.set_visible( False )
#
#        if i == sigma.shape[0]-1:
#            axs[i,j].set_xlabel( r'$r$ [km]' )
#        else:
#            axs[i,j].xaxis.set_visible( False )
#
#plt.subplots_adjust( hspace = 0.0, wspace = 0.0 )
#
##plt.savefig( 'fig.GaussianPerturbation.png', dpi = 300 )
#plt.show()
#plt.close()

# Make 2D plot

fig, axs = plt.subplots( sigma.shape[0], rho1.shape[0], figsize = (16,9) )

title  = r'$\rho(r)=\rho_0(r)+\rho_1\,e^{-(r-r_0)^2/(2\sigma^2)}\,\cos(\theta)$'
title += '\n'
title += r'$r_0=$ {:d} km'.format( np.int64( r0 ) )

fig.suptitle( title, fontsize = 15 )

xT, yT = 0.15, 0.7

ind = np.where( X1 > r0 )[0]
indR = ind[0]

Data = np.empty( (nDOFX[0],nDOFX[1]), np.float64 )

extent = [ xL[0], xU[0], xL[1], xU[1] ]

vmin = np.log10( 1.0e7 )
vmax = np.log10( 3.0e7 )
c = 'b'

for i in range( sigma.shape[0] ):

    for j in range( rho1.shape[0] ):

        rho11 = rho1[j] * rho0[indR]

        drho = rho11 * np.exp( -( X1 - r0 )**2 / ( 2.0 * sigma[i]**2 ) )

        for ii in range( nDOFX[1] ):

            Data[:,ii] = rho0 + drho * np.cos( X2[ii] )

        im = axs[i,j].imshow( np.log10( Data.T ), \
                              extent = extent, \
                              origin = 'lower', \
                              vmin = vmin, vmax = vmax, \
                              aspect = 'auto' )

        text = r'$\sigma$ = {:.0f} km'.format( sigma[i] )
        text += '\n'
        text += r'$\rho_1$ = {:.2f} $\rho_0(r_0)$'.format( rho1[j] )
        axs[i,j].text( xT, yT, text, transform = axs[i,j].transAxes, \
                       color = c )

        if j == 0:
            axs[i,j].set_ylabel( r'$\theta$' )
        else:
            axs[i,j].yaxis.set_visible( False )

        if i == sigma.shape[0]-1:
            axs[i,j].set_xlabel( r'$r$ [km]' )
        else:
            axs[i,j].xaxis.set_visible( False )

cax = fig.add_axes( [0.92,0.1,0.03,0.8] )
cbar = fig.colorbar( im, cax = cax )


plt.subplots_adjust( hspace = 0.0, wspace = 0.0 )

plt.savefig( 'fig.GaussianPerturbation.png', dpi = 300 )
#plt.show()
plt.close()

