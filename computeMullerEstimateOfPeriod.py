#!/usr/bin/env python3

import numpy as np

from TimeScales import TimeScales

rootDirectory = '/lump/data/accretionShockStudy/'

M    = np.array( [ '1.4', '2.0', '2.8' ], str )
Mdot = np.array( [ '0.3' ], str )
Rs   = np.array( [ '120', '150', '180' ], str )

TS = TimeScales()
T_Mul = np.empty( (M.shape[0],Rs.shape[0]), np.float64 )

for m in range( M.shape[0] ):
    for mdot in range( Mdot.shape[0] ):
        for rs in range( Rs.shape[0] ):

            ID \
              = 'GR1D_M{:}_Mdot{:}_Rs{:}'.format( M[m], Mdot[mdot], Rs[rs] )
            DataDirectory = rootDirectory + '{:}/'.format( ID )
            DataDirectory += '{:}.plt_00000000/'.format( ID )

            rInner = 4.0e1
            rOuter = np.float64( Rs[rs] )

            T_Mul[m,rs] \
              = TS.ComputeTimeScales( DataDirectory, rInner, rOuter )

np.savetxt( 'MulEst.dat', T_Mul )
