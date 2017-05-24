"""
Filenames of the tests
"""
import os
from os.path import join as pjoin

# This is a 4x4x4 q-point grid. The weights can be obtained from abinit.
nqpt = 8
wtq = [ 0.015625,  0.125   ,  0.0625  ,  0.09375 ,  0.375   ,  0.1875  ,
        0.046875,  0.09375 ]

dirname = os.path.dirname(__file__)

EIG0_fname = pjoin(dirname, 'EPC/WFK/out_data/odat_EIG.nc')

DDB_fnames = [
    pjoin(
        dirname,
        'DVSCF/qpt-{:0=4}/DVSCF/out_data/odat_DDB.nc'.format(iqpt+1)
        ) for iqpt in range(nqpt)
    ]

EIG_fnames = [
    pjoin(
        dirname,
        'EPC/qpt-{:0=4}/WFQ/out_data/odat_EIG.nc'.format(iqpt+1)
        ) for iqpt in range(nqpt)
    ]

EIGR2D_fnames = [
    pjoin(
        dirname,
        'EPC/qpt-{:0=4}/EPC/out_data/odat_EIGR2D.nc'.format(iqpt+1)
        ) for iqpt in range(nqpt)
    ]

GKK_fnames = [
    pjoin(
        dirname,
        'EPC/qpt-{:0=4}/EPC/out_data/odat_GKK.nc'.format(iqpt+1)
        ) for iqpt in range(nqpt)
    ]

fnames = dict(
    eigk_fname=EIG0_fname,
    eigq_fnames=EIG_fnames,
    ddb_fnames=DDB_fnames,
    eigr2d_fnames=EIGR2D_fnames,
    gkk_fnames=GKK_fnames,
    )

refdir = pjoin(dirname, 'epc_outputs')
