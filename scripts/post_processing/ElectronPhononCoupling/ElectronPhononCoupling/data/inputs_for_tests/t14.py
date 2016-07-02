from __future__ import print_function
import os
from os.path import join as pjoin

from ElectronPhononCoupling import compute_epc

# =========================================================================== #

DDB_fnames = """
../data_LiF/odat_calc_DS5_DDB.nc
../data_LiF/odat_calc_DS9_DDB.nc
../data_LiF/odat_calc_DS13_DDB.nc
""".split()

EIG_fnames = """
../data_LiF/odat_calc_DS6_EIG.nc
../data_LiF/odat_calc_DS10_EIG.nc
../data_LiF/odat_calc_DS14_EIG.nc
""".split()

EIGR2D_fnames = """
../data_LiF/odat_calc_DS7_EIGR2D.nc
../data_LiF/odat_calc_DS11_EIGR2D.nc
../data_LiF/odat_calc_DS15_EIGR2D.nc
""".split()

EIGI2D_fnames = """
../data_LiF/odat_calc_DS7_EIGI2D.nc
../data_LiF/odat_calc_DS11_EIGI2D.nc
../data_LiF/odat_calc_DS15_EIGI2D.nc
""".split()

GKK_fnames = """
../data_LiF/odat_calc_DS7_GKK.nc
../data_LiF/odat_calc_DS11_GKK.nc
../data_LiF/odat_calc_DS15_GKK.nc
""".split()

EIG0_fname = '../data_LiF/odat_calc_DS3_EIG.nc'


fnames = dict(
        eig0_fname=EIG0_fname,
        eigq_fnames=EIG_fnames,
        DDB_fnames=DDB_fnames,
        EIGR2D_fnames=EIGR2D_fnames,
        EIGI2D_fnames=EIGI2D_fnames,
        GKK_fnames=GKK_fnames,
        )

# =========================================================================== #

# This is a 2x2x2 q-point grid. The weights can be obtained from abinit.
nqpt = 3
wtq = [0.125, 0.5, 0.375]

epc = compute_epc(
        calc_type=1,
        write=True,
        output='output/t14',
        smearing_eV=0.01,
        temperature=True,
        temp_range=[0,600,300],
        lifetime=True,
        nqpt=nqpt,
        wtq=wtq,
        **fnames)


