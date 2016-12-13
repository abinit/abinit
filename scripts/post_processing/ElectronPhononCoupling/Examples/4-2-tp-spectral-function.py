"""
Compute the self-energy and the spectral function at various temperatures.
(ieig2rf=5)
"""

from ElectronPhononCoupling import compute_epc


# Lists of files used
# ===================

DDB_fnames = """
Calculations/01-LiF-dynamical/odat_calc_DS5_DDB.nc
Calculations/01-LiF-dynamical/odat_calc_DS9_DDB.nc
Calculations/01-LiF-dynamical/odat_calc_DS13_DDB.nc
""".split()

eigq_fnames = """
Calculations/01-LiF-dynamical/odat_calc_DS6_EIG.nc
Calculations/01-LiF-dynamical/odat_calc_DS10_EIG.nc
Calculations/01-LiF-dynamical/odat_calc_DS14_EIG.nc
""".split()

EIGR2D_fnames = """
Calculations/01-LiF-dynamical/odat_calc_DS7_EIGR2D.nc
Calculations/01-LiF-dynamical/odat_calc_DS11_EIGR2D.nc
Calculations/01-LiF-dynamical/odat_calc_DS15_EIGR2D.nc
""".split()

GKK_fnames = """
Calculations/01-LiF-dynamical/odat_calc_DS7_GKK.nc
Calculations/01-LiF-dynamical/odat_calc_DS11_GKK.nc
Calculations/01-LiF-dynamical/odat_calc_DS15_GKK.nc
""".split()

eig0_fname = 'Calculations/01-LiF-dynamical/odat_calc_DS3_EIG.nc'


# Computation of the self-energy and spectral function
# ====================================================

epc = compute_epc(
    calc_type = 4,          # Perform dynamical AHC calculation
    temperature = True,     # Do compute temperature dependence
    lifetime = False,       # Do not compute lifetime

    write = True,           # Do write the results
    output = 'Out/4-2',     # Rootname for the output
    
    smearing_eV = 0.05,                # Imaginary parameter for broadening.
    temp_range = [0, 800, 300],       # Temperature range (min, max, step)
    omega_range = [-0.1, 0.1, 0.0005], # Frequency range in Ha (min, max, step)

    nqpt = 3,                   # Number of q-points (2x2x2 qpt grid)
    wtq = [0.125, 0.5, 0.375],  # Weights of the q-points.
                                # These can be obtained by running Abinit
                                # with the corresponding k-point grid.
    
    eig0_fname = eig0_fname,        # All the files needed for
    eigq_fnames = eigq_fnames,      # this calculation.
    DDB_fnames = DDB_fnames,        #
    EIGR2D_fnames = EIGR2D_fnames,  #
    GKK_fnames = GKK_fnames,        # Note that instead of GKK.nc files,
                                    # one can also use the FAN.nc files
                                    # with the FAN_fnames keyword argument.
    )

