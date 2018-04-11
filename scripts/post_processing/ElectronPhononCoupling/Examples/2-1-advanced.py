"""
Initialize the EpcAnalyzer object and call specific functions
before writing the netCDF file. Only two functions are being called
in this example, but the user may look up the list of functions available.
"""

from ElectronPhononCoupling import EpcAnalyzer


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


# Initialization and computation
# ==============================

epca = EpcAnalyzer(

    rootname = 'Out/2-1',   # Rootname for the output
    
    smearing_eV = 0.10,               # Imaginary parameter for broadening.
    omega_range = [-0.1, 0.1, 0.001], # Frequency range in Ha (min, max, step)
    temp_range = [0, 1000, 250],      # Temperature range (min, max, step)

    nqpt = 3,                   # Number of q-points (2x2x2 qpt grid)
    wtq = [0.125, 0.5, 0.375],  # Weights of the q-points.
                                # These can be obtained by running Abinit
                                # with the corresponding k-point grid.
    
    eigk_fname = eigk_fname,        # All the files needed for
    eigq_fnames = eigq_fnames,      # this calculation.
    ddb_fnames = ddb_fnames,        #
    gkk_fnames = gkk_fnames,        #
    eigr2d_fnames = eigr2d_fnames,  #
    )


# Compute the temperature-dependent dynamical self-energy.
epca.compute_td_self_energy()

# Compute the mode-by-mode decomposition of the zero-point renormalization.
epca.compute_dynamical_zp_renormalization_modes()

# Write the result in a netCDF file.
epca.write_netcdf()
