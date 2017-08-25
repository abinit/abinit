"""
Compute the zero-point broadening (ZPB)
using the dynamical AHC theory (ieig2rf=5).
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


# Computation of the ZPB
# ======================

epc = compute(
    renormalization=False, # Do not compute the eigenvalues renormalization
    broadening = True,       # Do compute broadening
    temperature = True,    # Compute at several temperatures

    write = True,           # Do write the results
    rootname = 'Out/1-3',   # Rootname for the output

    smearing_eV = 0.01,         # Imaginary parameter for broadening.

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

