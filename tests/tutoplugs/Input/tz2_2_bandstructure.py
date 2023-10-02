#!/usr/bin/env python
from abipy.abilab import abiopen
import abipy.data as abidata
import matplotlib.pyplot as plt
import numpy as np

flist = ["tz2_2o_DS12_GSR.nc", "tz2_2o_DS22_GSR.nc"]
flist_gap = ["tz2_2o_DS13_FATBANDS.nc", "tz2_2o_DS23_FATBANDS.nc"]


gspec = dict(width_ratios=[2, 1], height_ratios=[1, 1])
fig, ax = plt.subplot_mosaic([['full0', 'gap0'], ['full5', 'gap5']],
                             gridspec_kw=gspec, figsize=(12, 8),
                             constrained_layout=True)
ax_list = ['full0', 'full5']
ax_gap_list = ['gap0', 'gap5']
for j in range(2):
    full = abiopen(flist[j]).ebands
    full.set_fermie_to_vbm()
    full.plot(with_gaps=False, ylims=(-3.0, 3.0), ax=ax[ax_list[j]], show=False)

    gap = abiopen(flist_gap[j]).ebands
    gap.set_fermie_to_vbm()
    gap.plot(with_gaps=True, ylims=(-1.0, 1.0), ax=ax[ax_gap_list[j]], show=False)


ax['full5'].set_title('5 GPa')
ax['full0'].set_title('0 GPa')
ax['gap5'].set_xticks([0, 17, len(gap.kpoints)-1])
ax['gap5'].set_xticklabels(['0.5H', 'A', '0.5L'])
ax['gap0'].set_xticks([0, 17, len(gap.kpoints)-1])
ax['gap0'].set_xticklabels(['0.5H', 'A', '0.5L'])

plt.savefig('tz2_2_bandstructure.pdf')
plt.show()
