"""
Compute the zero-point renormalization (ZPR) using the static AHC theory
(ieig2rf=1).
"""

from ElectronPhononCoupling import compute_epc


# Lists of files used
# ===================

DDB_fnames = """
Calculations/02-LiF-static/odat_calc_DS5_DDB.nc
Calculations/02-LiF-static/odat_calc_DS9_DDB.nc
Calculations/02-LiF-static/odat_calc_DS13_DDB.nc
""".split()

eigq_fnames = """
Calculations/02-LiF-static/odat_calc_DS6_EIG.nc
Calculations/02-LiF-static/odat_calc_DS10_EIG.nc
Calculations/02-LiF-static/odat_calc_DS14_EIG.nc
""".split()

EIGR2D_fnames = """
Calculations/02-LiF-static/odat_calc_DS7_EIGR2D.nc
Calculations/02-LiF-static/odat_calc_DS11_EIGR2D.nc
Calculations/02-LiF-static/odat_calc_DS15_EIGR2D.nc
""".split()

eig0_fname = 'Calculations/02-LiF-static/odat_calc_DS3_EIG.nc'


# Computation of the ZPR
# ======================

epc = compute_epc(
    calc_type = 1,          # Perform static AHC calculation
    temperature = False,    # Do not compute temperature dependence
    lifetime = False,       # Do not compute lifetime

    write = True,           # Do write the results
    output = 'Out/3-1',     # Rootname for the output
    
    nqpt = 3,                   # Number of q-points (2x2x2 qpt grid)
    wtq = [0.125, 0.5, 0.375],  # Weights of the q-points.
                                # These can be obtained by running Abinit
                                # with the corresponding k-point grid.
    
    eig0_fname = eig0_fname,        # All the files needed for
    eigq_fnames = eigq_fnames,      # this calculation.
    DDB_fnames = DDB_fnames,        #
    EIGR2D_fnames = EIGR2D_fnames,  #
    )

