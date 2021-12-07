#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

Meter    = 1.0
Kilogram = 1.0
Second   = 1.0

Kilometer   = 1.0e3      * Meter
SolarMass   = 1.98892e30 * Kilogram
Millisecond = 1.0e-3     * Second

G = 6.673e-11
c = 2.99792458e8

r = np.linspace( 40, 1000, 10000 ) * Kilometer

M_PNS = 1.4 * SolarMass

Mdot = 0.3 * SolarMass / Second # Accretion rate onto PNS

dt = 300.0 * Millisecond # Simulation time

dM_PNS = Mdot * dt

def LapseFunction( M ):
    return np.abs( ( 1.0 - M * G / ( 2.0 * c**2 * r ) ) \
                   / ( 1.0 + M * G / ( 2.0 * c**2 * r) ) )

def ConformalFactor( M ):
    return np.abs( 1.0 + M * G / ( 2.0 * c**2 * r ) )

alpha  = LapseFunction( M_PNS )
dalpha = LapseFunction( M_PNS + dM_PNS )

psi  = ConformalFactor( M_PNS )
dpsi = ConformalFactor( M_PNS + dM_PNS )

plt.plot( r / Kilometer, \
            np.abs( psi - dpsi ) / psi      , label = r'$|d\psi|/\psi$'     )
plt.plot( r / Kilometer, \
            np.abs( alpha - dalpha ) / alpha, label = r'$|d\alpha|/\alpha$' )
plt.xlabel( 'r [km]' )
plt.legend()
plt.show()

import os
os.system( 'rm -rf __pycache__ ' )
