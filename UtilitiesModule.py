#!/usr/bin/env python3

import numpy as np

def OverwriteFile( FileName, ForceChoice = False, Overwrite = False ):

    from os.path import isfile

    if( isfile( FileName ) ):

        Overwrite = False

    else:

        Overwrite = True

    return Overwrite

    if ForceChoice: return Overwrite

    Overwrite = True

    if( isfile( FileName ) ):

        YN = input( 'File: "{:}" exists. overwrite? (Y/N): '.format \
               ( FileName ) )

        if( not YN == 'Y' ):

            print( 'Not overwriting file' )
            Overwrite = False

    return Overwrite

# End of OverwriteFile


def GetFileArray( DataDirectory, PlotFileBaseName, Verbose = False ):

    from os import listdir

    if Verbose: print( '\nCalling GetFileArray...\n' )

    # Get last plotfile in directory

    FileArray = np.sort( np.array( [ f for f in listdir( DataDirectory ) ] ) )

    FileList = []

    for iFile in range( FileArray.shape[0] ):

        sFile = FileArray[iFile]

        if( sFile[0:len(PlotFileBaseName)+1] == PlotFileBaseName + '_' \
              and sFile[len(PlotFileBaseName)+1].isdigit() ):
            FileList.append( sFile )

    FileArray = np.array( FileList )

    if FileArray.shape[0] <= 0:

        msg += '\n  No files found in DataDirectory:'
        msg += ' {:s}\n\nDouble check the path\n'.format( DataDirectory )
        msg += '\n  Exiting...\n'
        assert FileArray.shape[0] > 0, msg

    return FileArray

# End of GetFileArray


def ChoosePlotFile \
      ( DataDirectory, PlotFileBaseName, argv = ['a'], \
        Verbose = False ):

    if Verbose: print( '\nCalling ChoosePlotFile...\n' )

    if( len( argv ) == 1 ):

        # Get last plotfile in directory

        FileArray = GetFileArray( DataDirectory, PlotFileBaseName, Verbose )
        File = FileArray[-1]

    elif( len( argv ) == 2 ):

        if argv[1][0].isalpha():

            File = argv[1]

        else:

            File = PlotFileBaseName + '_{:}'.format( argv[1].zfill(8) )

        if Verbose: print( '  File: {:}'.format( File ) )

    else:

        msg = 'Invalid number of optional parameters: {:d}'.format \
                ( len( argv ) )
        msg += '\n  Exiting...\n'
        assert len( argv ) <= 2, msg

    # Remove "/" at end of filename, if present
    if ( File[-1] == '/' ): File = File[:-1]

    return File

# End of ChoosePlotFile


def GetData( DataDirectory, PlotFileBaseName, argv, Field, Verbose = False ):

    import yt
    import numpy as np

    # https://yt-project.org/doc/faq/index.html#how-can-i-change-yt-s-log-level
    yt.funcs.mylog.setLevel(40) # Suppress yt warnings

    File = ChoosePlotFile( DataDirectory, PlotFileBaseName, argv, \
                           Verbose = Verbose )

    ds = yt.load( '{:}'.format( DataDirectory + File ) )

    MaxLevel = ds.index.max_level
    Time     = ds.current_time.to_ndarray()
    nX       = ds.domain_dimensions
    xL       = ds.domain_left_edge
    xU       = ds.domain_right_edge

    """
    https://yt-project.org/doc/reference/api/
    yt.data_objects.construction_data_containers.html#yt.data_objects.
    construction_data_containers.YTCoveringGrid
    """
    CoveringGrid \
      = ds.covering_grid \
          ( level           = MaxLevel, \
            left_edge       = xL, \
            dims            = nX * 2**MaxLevel, \
            num_ghost_zones = nX[0] )

    ds.force_periodicity()

    nDimsX = 1
    if nX[1] > 1: nDimsX += 1
    if nX[2] > 1: nDimsX += 1

    # --- Get Mesh ---

    xL = xL.to_ndarray()
    xU = xU.to_ndarray()

    dX = ( xU - xL ) / np.float64( nX )

    X1  = np.linspace( xL[0] + dX[0] / 2.0, xU[0] - dX[0] / 2.0, nX[0] )
    X2  = np.linspace( xL[1] + dX[1] / 2.0, xU[1] - dX[1] / 2.0, nX[1] )
    X3  = np.linspace( xL[2] + dX[2] / 2.0, xU[2] - dX[2] / 2.0, nX[2] )

    if  ( Field == 'PF_D'  ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = 'g/cm**3'

    elif( Field == 'PF_V1' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = 'km/s'

    elif( Field == 'PF_V2' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = '1/s'

    elif( Field == 'PF_V3' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = '1/s'

    elif( Field == 'PF_E'  ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = 'erg/cm**3'

    elif( Field == 'CF_D'  ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = 'g/cm**3'

    elif( Field == 'CF_S1' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = 'g/cm**2/s'

    elif( Field == 'CF_S2' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = 'g/cm/s'

    elif( Field == 'CF_S3' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = 'g/cm/s'

    elif( Field == 'CF_E'  ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = 'erg/cm**3'

    elif( Field == 'AF_P'  ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = 'erg/cm**3'

    elif( Field == 'AF_Cs' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = 'km/s'

    elif( Field == 'GF_Gm_11' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = ''

    elif( Field == 'GF_Gm_22' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = 'km**2'

    elif( Field == 'GF_Gm_33' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = 'km**2'

    elif( Field == 'GF_SqrtGm' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = 'km**2'

    elif( Field == 'GF_Psi' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = ''

    elif( Field == 'GF_Alpha' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = ''

    elif( Field == 'DF_TCI' ):

        Data = CoveringGrid[Field].to_ndarray()
        DataUnit = ''

    # --- Derived Fields ---

    elif( Field == 'pr4' ):

        p = CoveringGrid['AF_P'].to_ndarray()

        Data = np.empty( (nX[0],nX[1],nX[2]), np.float64 )

        for iX1 in range( nX[0] ):
            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):
                    Data[iX1,iX2,iX3] = p[iX1,iX2,iX3] * ( X1[iX1] * 1.0e5 )**4

        DataUnit = 'erg*cm'

    elif( Field == 'BernoulliConstant' ):

        c = 2.99792458e10

        rho   = CoveringGrid['PF_D'    ].to_ndarray()
        e     = CoveringGrid['PF_E'    ].to_ndarray()
        v1    = CoveringGrid['PF_V1'   ].to_ndarray() * 1.0e5
        v2    = CoveringGrid['PF_V2'   ].to_ndarray()
        p     = CoveringGrid['AF_P'    ].to_ndarray()
        alpha = CoveringGrid['GF_Alpha'].to_ndarray()
        Gm11  = CoveringGrid['GF_Gm11' ].to_ndarray()
        Gm22  = CoveringGrid['GF_Gm22' ].to_ndarray() * ( 1.0e5 )**2

        VSq = Gm11 * v1**2 + Gm22 * v2**2

        h = c**2 + ( e + p ) / rho
        W = 1.0 / np.sqrt( 1.0 - VSq / c**2 )

        B = alpha * h * W

        Data = B

        DataUnit = 'cm**2/s**2'

    elif( Field == 'BernoulliConstant_NR' ):

        c = 2.99792458e10

        rho   = CoveringGrid['PF_D'    ].to_ndarray()
        v1    = CoveringGrid['PF_V1'   ].to_ndarray() * 1.0e5
        v2    = CoveringGrid['PF_V2'   ].to_ndarray()
        p     = CoveringGrid['AF_P'    ].to_ndarray()
        Phi   = CoveringGrid['GF_Phi_N'].to_ndarray()
        Gm11  = CoveringGrid['GF_Gm_11' ].to_ndarray()
        Gm22  = CoveringGrid['GF_Gm_22' ].to_ndarray() * ( 1.0e5 )**2
        Gamma = CoveringGrid['AF_Gm'   ].to_ndarray()

        VSq = Gm11 * v1**2 + Gm22 * v2**2

        tau = Gamma / ( Gamma - 1.0 )

        B = 1.0 / 2.0 * VSq + tau * p / rho + Phi

        Data = B / c**2

        DataUnit = ''

    elif( Field == 'MassConstant_NR' ):

        rho = CoveringGrid['PF_D' ].to_ndarray()
        v   = CoveringGrid['PF_V1'].to_ndarray() * 1.0e5

        Data = X1**2 * rho * v

        DataUnit = 'g/s'

    elif( Field == 'PolytropicConstant' or Field == 'Entropy' ):

        PF_D  = CoveringGrid['PF_D' ].to_ndarray()
        AF_P  = CoveringGrid['AF_P' ].to_ndarray()
        AF_Gm = CoveringGrid['AF_Gm'].to_ndarray()

        Data  = AF_P / PF_D**AF_Gm

        DataUnit = 'erg/cm**3/(g/cm**3)**(4/3)'

    elif( Field == 'NonRelativisticSpecificEnthalpy' ):

        e   = CoveringGrid['PF_D'].to_ndarray()
        p   = CoveringGrid['PF_E'].to_ndarray()
        rho = CoveringGrid['AF_P'].to_ndarray()

        Data = ( e + p ) / rho

        DataUnit = 'cm**2/s**2'

    elif( Field == 'SpecificEnthalpy' ):

        c = 2.99792458e10

        e   = CoveringGrid['PF_D'].to_ndarray()
        p   = CoveringGrid['PF_E'].to_ndarray()
        rho = CoveringGrid['AF_P'].to_ndarray()

        Data = ( c**2 + ( e + p ) / rho ) / c**2

        DataUnit = ''

    elif( Field == 'LorentzFactor' ):

        c = 2.99792458e5

        Gm11 = CoveringGrid['GF_Gm_11'].to_ndarray()
        Gm22 = CoveringGrid['GF_Gm_22'].to_ndarray()
        Gm33 = CoveringGrid['GF_Gm_33'].to_ndarray()

        V1 = CoveringGrid['PF_V1'].to_ndarray()
        V2 = CoveringGrid['PF_V2'].to_ndarray()
        V3 = CoveringGrid['PF_V3'].to_ndarray()

        VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

        Data = 1.0 / np.sqrt( 1.0 - VSq / c**2 )

        DataUnit = ''

    elif( Field == 'TurbulentVelocity' ):

        Psi  = CoveringGrid['GF_Psi'  ].to_ndarray()
        Gm11 = CoveringGrid['GF_Gm_11'].to_ndarray()
        Gm22 = CoveringGrid['GF_Gm_22'].to_ndarray()
        Gm33 = CoveringGrid['GF_Gm_33'].to_ndarray()

        rho = CoveringGrid['PF_D' ].to_ndarray()
        V1  = CoveringGrid['PF_V1'].to_ndarray()
        V2  = CoveringGrid['PF_V2'].to_ndarray()
        V3  = CoveringGrid['PF_V3'].to_ndarray()

        # --- Compute angle-averaged and
        #     mass density weighted radial velocity ---

        AngleAveragedMass           = np.zeros( (nX[0]), np.float64 )
        AngleAveragedRadialVelocity = np.zeros( (nX[0]), np.float64 )

        Data = np.empty( nX, np.float64 )

        for iX1 in range( nX[0] ):

            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):

                    AngleAveragedMass[iX1] \
                      += rho[iX1,iX2,iX3] \
                           * Psi[iX1,iX2,iX3]**4 \
                           * np.sin( X2[iX2] ) * dX[1] * dX[2]

                    AngleAveragedRadialVelocity[iX1] \
                      += V1[iX1,iX2,iX3] * rho[iX1,iX2,iX3] \
                           * Psi[iX1,iX2,iX3]**4 \
                           * np.sin( X2[iX2] ) * dX[1] * dX[2]

            AngleAveragedRadialVelocity[iX1] /= AngleAveragedMass[iX1]

            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):

                    Data[iX1,iX2,iX3] \
                      = np.sqrt( \
                          Gm11[iX1,iX2,iX3] \
                            * ( V1[iX1,iX2,iX3] \
                                  - AngleAveragedRadialVelocity[iX1] )**2 \
                            + Gm22[iX1,iX2,iX3] * V2[iX1,iX2,iX3]**2 \
                            + Gm33[iX1,iX2,iX3] * V3[iX1,iX2,iX3]**2 )

        DataUnit = 'km/s'

    elif( Field == 'TurbulentEnergyDensity' ):

        Psi  = CoveringGrid['GF_Psi'  ].to_ndarray()
        Gm11 = CoveringGrid['GF_Gm_11'].to_ndarray()
        Gm22 = CoveringGrid['GF_Gm_22'].to_ndarray()
        Gm33 = CoveringGrid['GF_Gm_33'].to_ndarray()

        rho = CoveringGrid['PF_D' ].to_ndarray()
        V1  = CoveringGrid['PF_V1'].to_ndarray()
        V2  = CoveringGrid['PF_V2'].to_ndarray()
        V3  = CoveringGrid['PF_V3'].to_ndarray()

        AngleAveragedMass           = np.zeros( (nX[0]), np.float64 )
        AngleAveragedRadialVelocity = np.zeros( (nX[0]), np.float64 )

        c = 2.99792458e5

        Data = np.empty( nX, np.float64 )

        for iX1 in range( nX[0] ):

            # --- Compute angle-averaged and
            #     mass density weighted radial velocity ---

            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):

                    AngleAveragedMass[iX1] \
                      += rho[iX1,iX2,iX3] \
                           * Psi[iX1,iX2,iX3]**4 \
                           * np.sin( X2[iX2] ) * dX[1] * dX[2]

                    AngleAveragedRadialVelocity[iX1] \
                      += V1[iX1,iX2,iX3] * rho[iX1,iX2,iX3] \
                           * Psi[iX1,iX2,iX3]**4 \
                           * np.sin( X2[iX2] ) * dX[1] * dX[2]

            AngleAveragedRadialVelocity[iX1] /= AngleAveragedMass[iX1]

            # --- Compute turbulent energy density ---

            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):

                    # --- BetaSq = v_i * v^i / c^2 ---

                    BetaSq = ( Gm11[iX1,iX2,iX3] \
                                 * ( V1[iX1,iX2,iX3] \
                                       - AngleAveragedRadialVelocity[iX1] )**2 \
                                 + Gm22[iX1,iX2,iX3] * V2[iX1,iX2,iX3]**2 \
                                 + Gm33[iX1,iX2,iX3] * V3[iX1,iX2,iX3]**2 ) \
                               / c**2

                    W = 1.0 / np.sqrt( 1.0 - BetaSq )

                    Data[iX1,iX2,iX3] \
                      = rho[iX1,iX2,iX3] * ( c * 1.0e5 )**2 \
                          * W**2 * BetaSq / ( W + 1.0 )

        DataUnit = 'erg/cm**3'

    elif( Field == 'DivV2' ):

        # --- Sheck et al., (2008), A&A, 477, 931 ---

        V2 = CoveringGrid['PF_V2'].to_ndarray()

        Data = np.zeros( nX, np.float64 )

        indX2 = np.linspace( 1, nX[1]-2, nX[1]-2, dtype = np.int64 )

        for iX1 in range( nX[0] ):
            for iX2 in indX2:
                for iX3 in range( nX[2] ):

                    Data[iX1,iX2,iX3] \
                      = 1.0 / ( dX[1] * np.sin( X2[iX2] ) ) \
                          * (   np.sin( X2[iX2+1] ) * V2[iX1,iX2+1,iX3] \
                              - np.sin( X2[iX2-1] ) * V2[iX1,iX2-1,iX3] )

        for iX1 in range( nX[0] ):
            for iX3 in range( nX[2] ):

                Data[iX1,0 ,iX3] = Data[iX1,1 ,iX3]
                Data[iX1,-1,iX3] = Data[iX1,-2,iX3]


        DataUnit = 'rad/s'

    elif( Field == 'Vorticity' ):

        h1 = CoveringGrid['GF_h_1'].to_ndarray()
        h2 = CoveringGrid['GF_h_2'].to_ndarray()
        V1 = CoveringGrid['PF_V1' ].to_ndarray()
        V2 = CoveringGrid['PF_V2' ].to_ndarray()

        h1A = np.empty( (nX[0],nX[1]+2,nX[2]), np.float64 )
        h2A = np.empty( (nX[0],nX[1]+2,nX[2]), np.float64 )
        V1A = np.empty( (nX[0],nX[1]+2,nX[2]), np.float64 )
        V2A = np.empty( (nX[0],nX[1]+2,nX[2]), np.float64 )

        h1A[:,1:-1,:] = h1
        h2A[:,1:-1,:] = h2
        V1A[:,1:-1,:] = V1
        V2A[:,1:-1,:] = V2

        k = 0

        # --- Apply reflecting boundary conditions in theta ---

        for i in range( nX[0] ):

          h1A[i,0,k] = +h1A[i,1,k]
          h2A[i,0,k] = +h2A[i,1,k]
          V1A[i,0,k] = +V1A[i,1,k]
          V2A[i,0,k] = -V2A[i,1,k]

          h1A[i,nX[1]+1,k] = +h1A[i,nX[1],k]
          h2A[i,nX[1]+1,k] = +h2A[i,nX[1],k]
          V1A[i,nX[1]+1,k] = +V1A[i,nX[1],k]
          V2A[i,nX[1]+1,k] = -V2A[i,nX[1],k]

        # --- Compute vorticity in domain using
        #     central differences for derivatives ---

        Data = np.zeros( nX, np.float64 )

        k = 0
        for i in range( 1, nX[0] - 1 ):
            for j in range( 1, nX[1] + 1 ):

                Data[i,j-1,k] \
                  = 1.0 / ( h1A[i,j,k]      * h2A[i,j,k] ) \
                    * ( (   h2A[i+1,j,k]**2 * V2A[i+1,j,k] \
                          - h2A[i-1,j,k]**2 * V2A[i-1,j,k] ) \
                          / ( 2.0 * X1[i] ) \
                      - (   h1A[i,j+1,k]**2 * V1A[i,j+1,k] \
                          - h1A[i,j-1,k]**2 * V1A[i,j-1,k] ) \
                          / ( 2.0 * X2[j-1] ) )

        DataUnit = '1/s'

    else:

        print( '\nInvalid field: {:}'.format( Field ) )
        print( '\nValid choices:' )
        print( '--------------' )
        print( '  PF_D' )
        print( '  PF_V1' )
        print( '  PF_V2' )
        print( '  PF_V3' )
        print( '  PF_E' )
        print( '  CF_D' )
        print( '  CF_S1' )
        print( '  CF_S2' )
        print( '  CF_S3' )
        print( '  CF_E' )
        print( '  AF_P' )
        print( '  AF_Cs' )
        print( '  GF_Gm_11' )
        print( '  GF_Gm_22' )
        print( '  GF_Gm_33' )
        print( '  GF_Psi' )
        print( '  GF_Alpha' )
        print( '  DF_TCI' )
        print( '  pr4' )
        print( '  BernoulliConstant' )
        print( '  BernoulliConstant_NR' )
        print( '  MassConstant_NR' )
        print( '  PolytropicConstant' )
        print( '  NonRelativisticSpecificEnthalpy' )
        print( '  SpecificEnthalpy' )
        print( '  LorentzFactor' )
        print( '  TurbulentVelocity' )
        print( '  TurbulentEnergyDensity' )
        print( '  DivV2' )
        print( '  Vorticity' )
        exit( '\nExiting...' )

    if( nDimsX > 1 ):

        theta, r = np.meshgrid( X2, X1 )

        Data = Data[:,:,0]

    else:

        theta = X2
        r = X1
        Data = Data[:,0,0]

    return Data, DataUnit, r, theta, Time, xL, xU

# End of GetData


def GetNorm( UseLogScale, Data, vmin = +1.0e100, vmax = -1.0e100, \
             linthresh = 1.0e-2 ):

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, SymLogNorm

    if vmin > +1.0e99: vmin = Data.min()
    if vmax < -1.0e99: vmax = Data.max()

    if( UseLogScale ):

        if np.any( Data <= 0.0 ):

            Norm = SymLogNorm( vmin = vmin, vmax = vmax, \
                               linthresh = linthresh, base = 10 )

        else:

            Norm = LogNorm   ( vmin = vmin, vmax = vmax )

    else:

        Norm = plt.Normalize ( vmin = vmin, vmax = vmax )

    return Norm

# End of GetNorm
