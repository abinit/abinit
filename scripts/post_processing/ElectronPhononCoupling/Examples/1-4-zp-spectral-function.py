"""
Compute the self-energy and the spectral function at zero temperature.
(ieig2rf=5)
"""

from ElectronPhononCoupling import compute


# Lists of files used
# ===================

ddb_fnames = """
Calculations/01-LiF-dynamical/odat_calc_DS5_DDB.nc
Calculations/01-LiF-dynamical/odat_calc_DS9_DDB.nc
Calculations/01-LiF-dynamical/odat_calc_DS13_DDB.nc
""".split()

eigq_fnames = """
Calculations/01-LiF-dynamical/odat_calc_DS6_EIG.nc
Calculations/01-LiF-dynamical/odat_calc_DS10_EIG.nc
Calculations/01-LiF-dynamical/odat_calc_DS14_EIG.nc
""".split()

eigr2d_fnames = """
Calculations/01-LiF-dynamical/odat_calc_DS7_EIGR2D.nc
Calculations/01-LiF-dynamical/odat_calc_DS11_EIGR2D.nc
Calculations/01-LiF-dynamical/odat_calc_DS15_EIGR2D.nc
""".split()

gkk_fnames = """
Calculations/01-LiF-dynamical/odat_calc_DS7_GKK.nc
Calculations/01-LiF-dynamical/odat_calc_DS11_GKK.nc
Calculations/01-LiF-dynamical/odat_calc_DS15_GKK.nc
""".split()

eigk_fname = 'Calculations/01-LiF-dynamical/odat_calc_DS3_EIG.nc'


# Computation of the self-energy and spectral function
# ====================================================

epc = compute(

    renormalization=False,     # Do not compute the eigenvalues renormalization
    broadening = False,          # Do compute broadening
    self_energy = True,        # Compute frequency-dep. self-energy
    spectral_function = True,  # Compute frequency-dep. self-energy
    temperature = False,       # Compute only at T=0

    write = True,           # Do write the results
    rootname = 'Out/1-4',   # Rootname for the output
    
    smearing_eV = 0.05,               # Imaginary parameter for broadening.
    omega_range = [-0.1, 0.1, 0.001], # Frequency range in Ha (min, max, step)

    nqpt = 3,                   # Number of q-points (2x2x2 qpt grid)
    wtq = [0.125, 0.5, 0.375],  # Weights of the q-points.
                                # These can be obtained by running Abinit
                                # with the corresponding k-point grid.
    
    eigk_fname = eigk_fname,        # All the files needed for
    eigq_fnames = eigq_fnames,      # this calculation.
    ddb_fnames = ddb_fnames,        #
    eigr2d_fnames = eigr2d_fnames,  #
    gkk_fnames = gkk_fnames,        #
    )


# Plotting functions
# ==================

import netCDF4 as nc
import matplotlib.pyplot as plt

def plot_zp_self_energy(fname, ikpt=0, ieig=0):
    """Plot the self-energy computed at T=0."""

    with nc.Dataset(fname, 'r') as ds:
    
        # nomega
        omega = ds.variables['omegase'][...]
    
        # nspin, nkpt, nband, nomega
        self_energy = ds.variables['self_energy'][...]
    
    # Select one state
    re = self_energy[0,ikpt,ieig,:,0]
    im = self_energy[0,ikpt,ieig,:,1]
    
    plt.plot(omega, re, label='Real', linestyle='-')
    plt.plot(omega, im, label='Imag', linestyle='--')
    plt.legend(loc='upper right')

    fig = plt.gcf()
    plt.close()

    return fig


def plot_zp_spectral_function(fname, ikpt=0, ieig=0):
    """Plot the spectral function computed at T=0."""

    with nc.Dataset(fname, 'r') as ds:
    
        # nomega
        omega = ds.variables['omegase'][...]
    
        # nspin, nkpt, nband, nomega
        spectral_function = ds.variables['spectral_function'][...]
    
    # Select one state
    sf = spectral_function[0,ikpt,ieig,:]
    
    plt.plot(omega, sf)
    plt.ylim(0, 10)

    fig = plt.gcf()
    plt.close()

    return fig


# Do the plotting
# ================

fig_se = plot_zp_self_energy(epc.nc_output, ikpt=0, ieig=3)
fig_se.savefig('Out/1-4-plot-self-energy.eps')

fig_sf = plot_zp_spectral_function(epc.nc_output, ikpt=0, ieig=3)
fig_sf.savefig('Out/1-4-plot-spectral-function.eps')


