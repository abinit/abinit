"""
Filenames of the tests
"""
import os
from os.path import join as pjoin

from .. import LiF_g4


# This is a 2x2x2 q-point grid. The weights can be obtained from abinit.
nqpt = 3
wtq = [0.125, 0.5, 0.375]

# Indices of the q-points in the 4x4x4 grid.
iqpt_subset = [0, 2, 6]


dirname = os.path.dirname(__file__)

fnames = dict(
    eigk_fname=LiF_g4.fnames['eigk_fname'],
    eigq_fnames=list(),
    ddb_fnames=list(),
    eigr2d_fnames=list(),
    gkk_fnames=list(),
    )

for key in ('eigq_fnames', 'ddb_fnames', 'eigr2d_fnames', 'gkk_fnames'):
    for i in iqpt_subset:
        fnames[key].append(LiF_g4.fnames[key][i])

refdir = pjoin(dirname, 'epc_outputs')
