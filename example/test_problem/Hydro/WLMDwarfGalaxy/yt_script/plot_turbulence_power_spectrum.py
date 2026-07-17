import yt
import numpy as np
import matplotlib.pyplot as plt


def energy_spectrum_analytical(k):
    p_E     = -5./3.
    E_tot   =  2.0
    k_lower =  np.pi
    k_upper =  256*np.pi

    c_E = E_tot / ( 1./(p_E+1.) * ( (k_upper)**(p_E+1.) - (k_lower)**(p_E+1.) ) )
    return c_E * k**p_E


def get_scalar_power_spectrum(f, L_x1D, N_x1D, DCNormalized=False, Dimensionless=False, Total=False, VERBOSE=True):
    '''
    Parameters:
        f     : The scalar field
        L_x1D : Box size along each direction
        N_x1D : Number of cells along each direction
    '''

    V_tot = np.prod(L_x1D)                             # volume of the box
    N_tot = np.prod(N_x1D)                             # total number of cells = total number of k-modes

    N_selfc = np.prod(2-N_x1D%2)                       # number of self-conjugate (real) k-modes, including DC (and Nyquist when even)
    N_rednt = (N_tot-N_selfc)//2                       # number of redundant k-modes
    N_uqmin = (N_tot+N_selfc)//2                       # number of minimum unique k-modes
    N_rfftn = np.prod(N_x1D[:-1]) * (N_x1D[-1]//2 +1)  # number of k-modes output by real fft

    dx1D =  L_x1D / N_x1D                              # space interval along each direction
    dk1D =  2.0*np.pi / L_x1D                          # k-space interval along each direction

    k1D  = np.array( [2.0*np.pi*np.fft.fftfreq(N_x1D[d], d=dx1D[d]) for d in range(len(N_x1D))] )
    k1D_abs     = np.abs(k1D)                          # array of the magnitudes of the 1D k-modes
    k1D_min     = k1D_abs[k1D_abs>0].min()             # minimum magnitude of the 1D k-modes
    k1D_max     = k1D_abs[k1D_abs>0].max()             # maximum magnitude of the 1D k-modes
    n1D_Nyquist = N_x1D//2                             # index of the Nyquist mode

    # prepare the 3D k-space grid
    k3D_kx, k3D_ky, k3D_kz = np.meshgrid(k1D[0], k1D[1], k1D[2], indexing='ij')
    k3D_mag = np.sqrt(k3D_kx**2 + k3D_ky**2 + k3D_kz**2) # array of the magnitudes of the 3D k-modes
    k3D_min = k3D_mag[k3D_mag>0].min()                   # minimum magnitude of the 3D k-modes
    k3D_max = k3D_mag[k3D_mag>0].max()                   # maximum magnitude of the 3D k-modes

    # construct the field in x-space
    f_x3D        = f + 0j*f                              # complex field in x-space
    f_x3D_abs    = np.abs(f_x3D)                         # |f| field
    f_x3D_sq     = f_x3D_abs**2                          # |f|^2 field
    f_x3D_sq_sum = np.sum(f_x3D_sq)                      # sum of |f|^2
    f_x3D_avg    = np.average(f_x3D)                     # average of f
    f_x3D_min    = np.min(f_x3D_abs)                     # minimum of |f|
    f_x3D_max    = np.max(f_x3D_abs)                     # maximum of |f|

    # perform FFT to obtain the field in k-space
    f_k3D        = np.fft.fftn(f_x3D, norm='forward')    # f_k field in k-space, scaled by 1/N_tot in the forward transform
    f_k3D_abs    = np.abs(f_k3D)                         # |f_k| field
    f_k3D_sq     = f_k3D_abs**2                          # |f_k|^2 field
    f_k3D_sq_sum = np.sum(f_k3D_sq)                      # sum of |f_k|^2
    f_k3D_min    = np.min(f_k3D_abs)                     # minimum of |f_k|
    f_k3D_max    = np.max(f_k3D_abs)                     # maximum of |f_k|

    # DC mode
    DC = np.sqrt(f_k3D_sq.flat[0])                       # f_{k=0}

    # compute the power
    P_k3D      = f_k3D_sq * V_tot                        # |f_k|^2 * V
    P_k3D_sum  = np.sum(P_k3D)                           # sum of |f_k|^2 * V
    P_k3D_min  = np.min(P_k3D)                           # maximum of |f_k|^2 * V
    P_k3D_max  = np.max(P_k3D)                           # minimum of |f_k|^2 * V

    # compute the second moment of the field in x-space and k-space
    SecondMoment_x3D = f_x3D_sq_sum / N_tot              # expectation value of |f|^2
    SecondMoment_k3D = P_k3D_sum    / V_tot              # expectation value of |f|^2

    # prepare the bins for the power spectrum
    kBin_dk          = np.min(dk1D)
    kBin_center_min  = 0.0                                                                          # minimum of k-bins
    kBin_center_max  = k3D_max                                                                      # maximum of k-bins
    kBin_edge        = np.arange(kBin_center_min-0.5*kBin_dk, kBin_center_max+1.0*kBin_dk, kBin_dk) # k-bin edges
    kBin_center      = kBin_edge[1:] - 0.5*kBin_dk                                                  # k-bin centers

    # calculate the power spectrum
    NumModes_kBin, _ = np.histogram(k3D_mag, kBin_edge)                                         # number of k-modes in each bin
    TotPower_kBin, _ = np.histogram(k3D_mag, kBin_edge, weights=P_k3D)                          # sum of P(k) in each bin
    AvgPower_kBin    = TotPower_kBin / NumModes_kBin                                            # 1D binned power spectrum P(k_b)
    IntPower_kBin    = 4.0 * np.pi * kBin_center**2 * AvgPower_kBin * kBin_dk / (2.0*np.pi)**3  # integral of power spectrum

    DCNormalizedPower_kBin  = AvgPower_kBin / DC**2                                              # scale the DC to 1, and the power spectrum becomes in units of volume
    DimensionlessPower_kBin = DCNormalizedPower_kBin * 4.0*np.pi*kBin_center**3/(2.0*np.pi)**3   # dimensionless form used in cosmological matter power spectrum

    if VERBOSE:
        np.set_printoptions(formatter={'float': '{: 12.4e}'.format})
        print('')
        print(f'Box size along x               = {L_x1D[0]            = : 12.4e}')
        print(f'Box size along y               = {L_x1D[1]            = : 12.4e}')
        print(f'Box size along z               = {L_x1D[2]            = : 12.4e}')
        print(f'Box volume                     = {V_tot               = : 12.4e}')

        print('')
        print(f'Number of cells x              = {N_x1D[0]            = : 12d}')
        print(f'Number of cells y              = {N_x1D[1]            = : 12d}')
        print(f'Number of cells z              = {N_x1D[2]            = : 12d}')
        print(f'Total number of cells          = {N_tot               = : 12d}')

        print('')
        print(f'Number of self-conjugate modes = {N_selfc             = : 12d}')
        print(f'Number of redundant modes      = {N_rednt             = : 12d}')
        print(f'Number of minimum unique modes = {N_uqmin             = : 12d}')
        print(f'Number of real fft modes       = {N_rfftn             = : 12d}')

        print('')
        print(f'Interval in x-space            = {dx1D                = }')
        print(f'Interval in k-space            = {dk1D                = }')

        print('')
        print(f'{k1D.shape           = }')
        print(f'{k1D[0].shape        = }, {k1D[0][0] = : 12.4e}, {k1D[0][1] = : 12.4e}, {k1D[0][n1D_Nyquist[0]-1] = : 12.4e}, {k1D[0][n1D_Nyquist[0]] = : 12.4e}, {k1D[0][-1] = : 12.4e}')
        print(f'{k1D[1].shape        = }, {k1D[1][0] = : 12.4e}, {k1D[1][1] = : 12.4e}, {k1D[1][n1D_Nyquist[1]-1] = : 12.4e}, {k1D[1][n1D_Nyquist[1]] = : 12.4e}, {k1D[1][-1] = : 12.4e}')
        print(f'{k1D[2].shape        = }, {k1D[2][0] = : 12.4e}, {k1D[2][1] = : 12.4e}, {k1D[2][n1D_Nyquist[2]-1] = : 12.4e}, {k1D[2][n1D_Nyquist[2]] = : 12.4e}, {k1D[2][-1] = : 12.4e}')
        print(f'{k1D_min             = : 12.4e}')
        print(f'{k1D_max             = : 12.4e}')
        print(f'{n1D_Nyquist[0]      = : 12d}')
        print(f'{n1D_Nyquist[1]      = : 12d}')
        print(f'{n1D_Nyquist[2]      = : 12d}')

        print('')
        print(f'{k3D_kx.shape        = }')
        print(f'{k3D_ky.shape        = }')
        print(f'{k3D_kz.shape        = }')

        print('')
        print(f'{k3D_mag.shape       = }')
        print(f'{k3D_min             = : 12.4e}')
        print(f'{k3D_max             = : 12.4e}')

        print('')
        print(f'{f_x3D.shape         = }')
        print(f'{f_x3D_sq_sum        = : 12.4e}')
        print(f'{f_x3D_avg           = : 12.4e}')
        print(f'{f_x3D_min           = : 12.4e}')
        print(f'{f_x3D_max           = : 12.4e}')

        print('')
        print(f'{f_k3D.shape         = }')
        print(f'{f_k3D_sq_sum *N_tot = : 12.4e}')
        print(f'{f_k3D_min           = : 12.4e}')
        print(f'{f_k3D_max           = : 12.4e}')

        print('')
        print(f'{DC                  = : 12.4e}')

        print('')
        print(f'{P_k3D.shape         = }')
        print(f'{P_k3D_sum           = : 12.4e}')
        print(f'{P_k3D_min           = : 12.4e}')
        print(f'{P_k3D_max           = : 12.4e}')
        print('')
        print(f'{SecondMoment_x3D    = : 12.4e}')
        print(f'{SecondMoment_k3D    = : 12.4e}')

        print('')
        print(f'{kBin_dk             = : 12.4e}')
        print(f'{kBin_center_min     = : 12.4e}')
        print(f'{kBin_center_max     = : 12.4e}')
        print(f'{kBin_edge.shape     = }')
        print(f'{kBin_center.shape   = }')
        print(f'kBin_edge           = [{kBin_edge[ 0] : 12.4e}, {kBin_edge[ 1] : 12.4e}, {kBin_edge[ 2] : 12.4e},       ...,       {kBin_edge[-2] : 12.4e}, {kBin_edge[-1] : 12.4e}]')
        print(f'kBin_center         = [      {kBin_center[ 0] : 12.4e}, {kBin_center[ 1] : 12.4e}, {kBin_center[ 2] : 12.4e}, ..., {kBin_center[-2] : 12.4e}, {kBin_center[-1] : 12.4e}      ]')

        print('')
        print(f'{NumModes_kBin.shape = }')
        print(f'{TotPower_kBin.shape = }')
        print(f'{AvgPower_kBin.shape = }')
        print(f'NumModes_kBin       = [      {NumModes_kBin[ 0] : 12d}, {NumModes_kBin[ 1] : 12d}, {NumModes_kBin[ 2] : 12d}, ..., {NumModes_kBin[-2] : 12d}, {NumModes_kBin[-1] : 12d}      ]')
        print(f'TotPower_kBin       = [      {TotPower_kBin[ 0] : 12.4e}, {TotPower_kBin[ 1] : 12.4e}, {TotPower_kBin[ 2] : 12.4e}, ..., {TotPower_kBin[-2] : 12.4e}, {TotPower_kBin[-1] : 12.4e}      ]')
        print(f'AvgPower_kBin       = [      {AvgPower_kBin[ 0] : 12.4e}, {AvgPower_kBin[ 1] : 12.4e}, {AvgPower_kBin[ 2] : 12.4e}, ..., {AvgPower_kBin[-2] : 12.4e}, {AvgPower_kBin[-1] : 12.4e}      ]')
        print(f'IntPower_kBin       = [      {IntPower_kBin[ 0] : 12.4e}, {IntPower_kBin[ 1] : 12.4e}, {IntPower_kBin[ 2] : 12.4e}, ..., {IntPower_kBin[-2] : 12.4e}, {IntPower_kBin[-1] : 12.4e}      ]')
        print(f'{sum(NumModes_kBin)  = : 12d}')
        print(f'{sum(TotPower_kBin)  = : 12.4e}')
        print(f'{sum(IntPower_kBin)  = : 12.4e}') # this should reproduce the second moment

    # output the results
    if Dimensionless:
        return kBin_center, DimensionlessPower_kBin

    elif DCNormalized:
        return kBin_center, DCNormalizedPower_kBin

    elif Total:
        return kBin_center, TotPower_kBin

    else:
        return kBin_center, AvgPower_kBin


def get_vector_power_spectrum(v, L_x1D, N_x1D, Total=False):

    kBin_center_vx, Power_kBin_vx = get_scalar_power_spectrum(v[0], L_x1D, N_x1D, Total=Total)
    kBin_center_vy, Power_kBin_vy = get_scalar_power_spectrum(v[1], L_x1D, N_x1D, Total=Total)
    kBin_center_vz, Power_kBin_vz = get_scalar_power_spectrum(v[2], L_x1D, N_x1D, Total=Total)

    kBin_center = kBin_center_vx
    Power_kBin  = Power_kBin_vx + Power_kBin_vy + Power_kBin_vz

    return kBin_center, Power_kBin


def get_specific_kinetic_energy_power_spectrum(VelocityField, L_x1D, N_x1D):
    kBin_center, TotPower_kBin = get_vector_power_spectrum(VelocityField, L_x1D, N_x1D, Total=True)

    deltaK_kBin = np.gradient(kBin_center)
    V_tot = np.prod(L_x1D)

    SpecificKineticEnergyPower_kBin = 0.5*TotPower_kBin/(V_tot * deltaK_kBin)

    print(f'{np.sum(SpecificKineticEnergyPower_kBin*deltaK_kBin)=}')

    return kBin_center, SpecificKineticEnergyPower_kBin


def get_velocity_field(L_x1D, N_x1D, random_seed=123):
    """
    Sample a zero-mean, homogeneous, isotropic, Gaussian, incompressible 3D velocity field.
    Inputs: box sizes Lx,Ly,Lz; grid sizes Nx,Ny,Nz; analytic specific kinetic-energy spectrum E(k).
    Returns: real-space velocity grid with shape (3, Nz, Ny, Nx).
    """

    Lx     = L_x1D[0]
    Ly     = L_x1D[1]
    Lz     = L_x1D[2]

    Nx     = N_x1D[0]
    Ny     = N_x1D[1]
    Nz     = N_x1D[2]

    V_tot  = np.prod(L_x1D)        # total box volume
    N_tot  = np.prod(N_x1D)        # total number of grids

    # set up the k-space
    k3D_kx, k3D_ky, k3D_kz = np.meshgrid( 2 * np.pi * np.fft.fftfreq( Nx, Lx / Nx), # grids of the components of wavevector k
                                          2 * np.pi * np.fft.fftfreq( Ny, Ly / Ny), # includes positive and negative wavenumber along x and y
                                          2 * np.pi * np.fft.rfftfreq(Nz, Lz / Nz), # only positive wavenumber, using rfftn layout, along z
                                          indexing="ij" )                           # shape = (Nx, Ny, Nz)
    kVector_k3D            = np.array([k3D_kx, k3D_ky, k3D_kz])                     # k vector field; shape=(3, Nx, Ny, Nz)

    k_sq_k3D               = (kVector_k3D * kVector_k3D).sum(0)                     # |k|^2 = k_x^2 + k_y^2 + k_z^2; shape=(Nx, Ny, Nz)
    nonDC                  = k_sq_k3D > 0

    # draw a random 3-component complex Gaussian vector field
    # --> each compoent is (X+iY)/sqrt(2) with X,Y ~ N(mu=0, variance=1)
    # --> total variance is 3
    rng = np.random.default_rng(seed=random_seed)
    realGaussianVector_k3D    = rng.normal(size=kVector_k3D.shape)
    imagGaussianVector_k3D    = rng.normal(size=kVector_k3D.shape)
    complexGaussianVector_k3D = (realGaussianVector_k3D + 1j * imagGaussianVector_k3D) / np.sqrt(2)
    print(f'{np.mean((np.abs(complexGaussianVector_k3D)**2).sum(0)) = }')

    # impose the incompressible fluid condition
    # --> div(v) = 0 -> k dot \tilde{v} = 0
    # --> project the vector field onto the plane perpendicular to k
    # --> total variance becomes 2
    complexGaussianVector_kParallel_k3D      = kVector_k3D * (kVector_k3D * complexGaussianVector_k3D).sum(0, keepdims=True) / np.where(nonDC, k_sq_k3D, 1)
    complexGaussianVector_kPerpendicular_k3D = complexGaussianVector_k3D - complexGaussianVector_kParallel_k3D
    print(f'{np.mean((np.abs(complexGaussianVector_kParallel_k3D)**2).sum(0)) = }')
    print(f'{np.mean((np.abs(complexGaussianVector_kPerpendicular_k3D)**2).sum(0)) = }')

    # compute the amplitude to scale the variance at different k
    # --> DC mode is zero
    # --> non-DC mode follows the input specific kinetic energy spectrum
    amplitude_k3D        = np.zeros_like(k_sq_k3D)
    amplitude_k3D[nonDC] = np.sqrt(V_tot * 2 * np.pi**2 * energy_spectrum_analytical(np.sqrt(k_sq_k3D[nonDC])) / k_sq_k3D[nonDC])

    # velocity field in the Fourier space
    vVector_k3D = amplitude_k3D * complexGaussianVector_kPerpendicular_k3D

    # enforce Hermitian symmetry, F_{-n_x, -n_y, -n_z} = F_{n_x, n_y, n_z}^*
    # for most of n_z, the conjugate partner is not stored since we use irfftn()
    # we only need to deal with the planes where conjugate partners exist
    # --> i.e., n_z is self-conjugate: n_z = -n_z mod Nz
    # --> (1) n_z = 0, DC always, and (2) n_z = n_z_Nyquist when Nz is even
    n_z_selfConjugate = [0, -1] if Nz%2 == 0 else [0]
    for n_z in n_z_selfConjugate:
        # for each (n_x, n_y, n_z), average it with its Hermitian conjugate partner (-n_x, -n_y, n_z) to make it Hermitian symmetric
        # normalized by sqrt(2) to perserve the variance
        n_x_partner = ( -np.arange(Nx) ) % Nx
        n_y_partner = ( -np.arange(Ny) ) % Ny
        vVector_k3D[:, :, :, n_z] = (vVector_k3D[:, :, :, n_z] + np.conj(vVector_k3D[:, n_x_partner][:, :, n_y_partner, n_z])) / np.sqrt(2)

    # convert the continuous-space Fourier coefficients into DFT coefficients
    dV              = V_tot / N_tot     # cell volume
    vVector_DFT_k3D = vVector_k3D / dV  #

    # perform the inverse Fourier transform
    vVector_x3D = np.fft.irfftn(vVector_DFT_k3D, s=(Nx, Ny, Nz), axes=(1, 2, 3)).real  # real-space velocity vector field

    return vVector_x3D


def plot_power_spectrum(k_plot, Spectrum, k1D_Nyquist=None, PLOT_ANALYTICAL=True):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if k1D_Nyquist is not None:
        ax.axvline(k1D_Nyquist, color='grey')

    ax.loglog( k_plot, Spectrum, '-o', markersize=3, color='b', label='Power Spectrum' )

    if PLOT_ANALYTICAL:
        ax.loglog( k_plot, energy_spectrum_analytical(k_plot), '-.', color='r', label='Analytical' )
    else:
        ax.loglog( k_plot, Spectrum[5]*(k_plot/k_plot[5])**(-2), '-.', color='r', label='Analytical' )


    ax.set_xlabel( r'$k$' )
    ax.set_ylabel( r'$P(k)$' )
    ax.grid()
    ax.legend()

    plt.tight_layout()
    fig.savefig( 'fig_power_spectrum.png' )
    plt.close()


def main():

    L_x1D       = np.array([2.0, 2.0, 2.0])
    N_x1D       = np.array([128, 128, 128], dtype=int)
    k1D_Nyquist = np.pi / np.min(L_x1D/N_x1D)

    velocity_field = get_velocity_field(L_x1D, N_x1D)

    k_plot, Spectrum = get_specific_kinetic_energy_power_spectrum(velocity_field, L_x1D, N_x1D)

    nonDC = slice(1, None)
    plot_power_spectrum(k_plot[nonDC], Spectrum[nonDC], k1D_Nyquist=k1D_Nyquist)

    table = np.genfromtxt( 'Tur_Table.dat', delimiter=None, comments='#',
                           names=['x','y','z','vx','vy','vz'], dtype=None, encoding=None )
    velocity_field_table = np.array( [table['vx'].reshape(129,129,129), table['vy'].reshape(129,129,129), table['vz'].reshape(129,129,129)] )
    print(table['vx'].shape)
    print(table['vy'].shape)
    print(table['vz'].shape)
    print(velocity_field_table.shape)

    k_plot, Spectrum = get_specific_kinetic_energy_power_spectrum(velocity_field_table, np.array([2.0*np.pi, 2.0*np.pi, 2.0*np.pi]), np.array([129, 129, 129], dtype=int))

    nonDC = slice(1, None)
    plot_power_spectrum(k_plot[nonDC], Spectrum[nonDC], k1D_Nyquist=64.0, PLOT_ANALYTICAL=False)



if __name__ == '__main__':
    main()
