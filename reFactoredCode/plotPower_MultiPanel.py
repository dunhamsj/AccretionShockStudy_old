#!/usr/bin/env python3

import numpy as np
from os.path import isdir
import matplotlib.pyplot as plt
from FitPowerToModel import FittingFunction

R    = np.array( [ 'NR', 'GR' ], str )
M    = np.array( [ '1.4', '2.0', '2.8' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '120', '150', '180' ], str )

ID = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )

t  = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )
P0 = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )
P1 = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )
P2 = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )
P3 = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )
P4 = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )

t0      = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )
t1      = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )
LogF    = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )
omegaR  = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )
omegaI  = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )
delta   = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )
domegaR  = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )
domegaI  = np.empty( (R.shape[0],M.shape[0],Mdot.shape[0],Rs.shape[0]), object )

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

for m in range( M.shape[0] ):
    for rs in range( Rs.shape[0] ):
        axs[m,rs].semilogy( t[0,m,0,rs], P1[0,m,0,rs] )
        axs[m,rs].semilogy( t[1,m,0,rs], P1[1,m,0,rs] )

import os
os.system( 'rm -rf __pycache__ ' )
