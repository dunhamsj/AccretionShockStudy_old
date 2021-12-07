#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

class Tally:


    def __init__( self, DataDirectory, ID ):

        self.DataDirectory = DataDirectory
        self.ID            = ID

        return


    def GetData( self ):

        Root = self.DataDirectory + self.ID + '/' + self.ID \
                 + '.plt_StandingAccretionShock_Relativis_Tally_'

        self.Mass   = np.loadtxt( Root + 'BaryonicMass.dat', skiprows = 1 )
        self.Energy = np.loadtxt( Root + 'Energy.dat', skiprows = 1 )


    def PlotData( self ):

        self.GetData()

        fig, axs = plt.subplots( 2, 1 )

        axs[0].plot( self.Mass  [:,0], self.Mass  [:,4] )
        axs[1].plot( self.Energy[:,0], self.Energy[:,4] )

        axs[0].set_ylabel( 'Change [Msun]' )
        axs[1].set_ylabel( 'Change [B]' )

        plt.show()

if __name__ == '__main__':

    DataDirectory \
      = '/Users/dunhamsj/Research/Data/AccretionShockParameterStudy/'

    ID = 'GR2D_M1.4_Mdot0.3_Rs180_PA0.004_nX320x064_smooth'

    T = Tally( DataDirectory, ID )

    T.PlotData()

    import os
    os.system( 'rm -rf __pycache__ ' )
