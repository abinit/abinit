
import os
from os.path import join as pjoin

dirname = os.path.dirname(__file__)

nqpt = 3
wtq = [0.125, 0.5, 0.375]

DDB_fnames = [ pjoin(dirname, fname) for fname in """
                    odat_calc_DS5_DDB.nc
                    odat_calc_DS9_DDB.nc
                    odat_calc_DS13_DDB.nc
                    """.split() ]

EIG_fnames = [ pjoin(dirname, fname) for fname in """
                    odat_calc_DS6_EIG.nc
                    odat_calc_DS10_EIG.nc
                    odat_calc_DS14_EIG.nc
                    """.split() ]

EIGR2D_fnames = [ pjoin(dirname, fname) for fname in """
                    odat_calc_DS7_EIGR2D.nc
                    odat_calc_DS11_EIGR2D.nc
                    odat_calc_DS15_EIGR2D.nc
                   """.split() ]

EIGI2D_fnames = [ pjoin(dirname, fname) for fname in """
                    odat_calc_DS7_EIGI2D.nc
                    odat_calc_DS11_EIGI2D.nc
                    odat_calc_DS15_EIGI2D.nc
                   """.split() ]

FAN_fnames = [ pjoin(dirname, fname) for fname in """
                    odat_calc_DS7_FAN.nc
                    odat_calc_DS11_FAN.nc
                    odat_calc_DS15_FAN.nc
                   """.split() ]

GKK_fnames = [ pjoin(dirname, fname) for fname in """
                    odat_calc_DS7_GKK.nc
                    odat_calc_DS11_GKK.nc
                    odat_calc_DS15_GKK.nc
                   """.split() ]

EIG0_fname = pjoin(dirname, 'odat_calc_DS3_EIG.nc')


fnames = dict(
        eigk_fname=EIG0_fname,
        eigq_fnames=EIG_fnames,
        ddb_fnames=DDB_fnames,
        eigr2d_fnames=EIGR2D_fnames,
        eigi2d_fnames=EIGI2D_fnames,
        #fan_fnames=FAN_fnames,
        gkk_fnames=GKK_fnames,
        )

refdir = pjoin(dirname, 'epc_outputs')
