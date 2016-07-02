
import os
from os.path import join as pjoin

LiF_datadir = pjoin(os.path.dirname(__file__), 'data_LiF')

LiF_nqpt = 3
LiF_wtq = [0.125, 0.5, 0.375]

LiF_DDB_fnames = [ pjoin(LiF_datadir, fname) for fname in """
                    odat_calc_DS5_DDB.nc
                    odat_calc_DS9_DDB.nc
                    odat_calc_DS13_DDB.nc
                    """.split() ]

LiF_EIG_fnames = [ pjoin(LiF_datadir, fname) for fname in """
                    odat_calc_DS6_EIG.nc
                    odat_calc_DS10_EIG.nc
                    odat_calc_DS14_EIG.nc
                    """.split() ]

LiF_EIGR2D_fnames = [ pjoin(LiF_datadir, fname) for fname in """
                    odat_calc_DS7_EIGR2D.nc
                    odat_calc_DS11_EIGR2D.nc
                    odat_calc_DS15_EIGR2D.nc
                   """.split() ]

LiF_EIGI2D_fnames = [ pjoin(LiF_datadir, fname) for fname in """
                    odat_calc_DS7_EIGI2D.nc
                    odat_calc_DS11_EIGI2D.nc
                    odat_calc_DS15_EIGI2D.nc
                   """.split() ]

LiF_FAN_fnames = [ pjoin(LiF_datadir, fname) for fname in """
                    odat_calc_DS7_FAN.nc
                    odat_calc_DS11_FAN.nc
                    odat_calc_DS15_FAN.nc
                   """.split() ]

LiF_GKK_fnames = [ pjoin(LiF_datadir, fname) for fname in """
                    odat_calc_DS7_GKK.nc
                    odat_calc_DS11_GKK.nc
                    odat_calc_DS15_GKK.nc
                   """.split() ]

LiF_EIG0_fname = pjoin(LiF_datadir, 'odat_calc_DS3_EIG.nc')


LiF_fnames = dict(
        eig0_fname=LiF_EIG0_fname,
        eigq_fnames=LiF_EIG_fnames,
        DDB_fnames=LiF_DDB_fnames,
        EIGR2D_fnames=LiF_EIGR2D_fnames,
        EIGI2D_fnames=LiF_EIGI2D_fnames,
        #FAN_fnames=LiF_FAN_fnames,
        GKK_fnames=LiF_GKK_fnames,
        )

LiF_outputdir = pjoin(os.path.dirname(__file__), 'outputs_of_tests')

LiF_outputs = dict()

for t in ['t11', 't12', 't13', 't14',
          't21', 't22', 't23', 't24',
          't31', 't32', 't33', 't34']:
    LiF_outputs[t] = pjoin(LiF_outputdir, t + '_EP.nc')





