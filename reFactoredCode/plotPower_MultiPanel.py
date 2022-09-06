#!/usr/bin/env python3

import numpy as np
from os.path import isdir
import matplotlib.pyplot as plt
from FitPowerToModel import FittingFunction

R    = np.array( [ 'NR', 'GR' ], str )
M    = np.array( [ '1.4', '2.0', '2.8' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '120', '150', '180' ], str )

arrShape = (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0])

ID = np.empty( arrShape, object )

t  = np.empty( arrShape, object )
P0 = np.empty( arrShape, object )
P1 = np.empty( arrShape, object )
P2 = np.empty( arrShape, object )
P3 = np.empty( arrShape, object )
P4 = np.empty( arrShape, object )

t0      = np.empty( arrShape, object )
t1      = np.empty( arrShape, object )
LogF    = np.empty( arrShape, object )
omegaR  = np.empty( arrShape, object )
omegaI  = np.empty( arrShape, object )
delta   = np.empty( arrShape, object )
domegaR = np.empty( arrShape, object )
domegaI = np.empty( arrShape, object )

for r in range( R.shape[0] ):
    for m in range( M.shape[0] ):
        for mdot in range( Mdot.shape[0] ):
            for rs in range( Rs.shape[0] ):

                ID[r,m,mdot,rs] \
                  = '{:}2D_M{:}_Mdot{:}_Rs{:}'.format \
                     ( R[r], M[m], Mdot[mdot], Rs[rs] )

                plotFileDirectory \
                  = '/lump/data/accretionShockStudy/{:}/'.format \
                    ( ID[r,m,mdot,rs] )

                if not isdir( plotFileDirectory ): continue

                plotFileBaseName = '{:}.plt_'.format( ID[r,m,mdot,rs] )

                dataFileName \
                  = '.{:}_LegendrePowerSpectrum.dat'.format( ID[r,m,mdot,rs] )

                t [r,m,mdot,rs], \
                P0[r,m,mdot,rs], \
                P1[r,m,mdot,rs], \
                P2[r,m,mdot,rs], \
                P3[r,m,mdot,rs], \
                P4[r,m,mdot,rs] \
                  = np.loadtxt( dataFileName )

                dataFileName = '.{:}_Fit.dat'.format( ID[r,m,mdot,rs] )

                t0        [r,m,mdot,rs], \
                t1        [r,m,mdot,rs], \
                LogF      [r,m,mdot,rs], \
                omegaR    [r,m,mdot,rs], \
                omegaI    [r,m,mdot,rs], \
                delta     [r,m,mdot,rs], \
                dummy0, \
                omegaR_err[r,m,mdot,rs], \
                omegaI_err[r,m,mdot,rs], \
                dummy1 \
                  = np.loadtxt( dataFileName )

fig, axs = plt.subplots( M.shape[0], Rs.shape[0], figsize = (12,9) )

nr   = 0
gr   = 1
mdot = 0
for m in range( M.shape[0] ):
    for rs in range( Rs.shape[0] ):
        axs[m,rs].semilogy( t[nr,m,mdot,rs], P1[nr,m,mdot,rs] )
        axs[m,rs].semilogy( t[gr,m,mdot,rs], P1[gr,m,mdot,rs] )

import os
os.system( 'rm -rf __pycache__ ' )
