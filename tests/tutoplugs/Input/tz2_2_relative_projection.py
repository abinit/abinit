#!/usr/bin/env python
from abipy.abilab import abiopen
import abipy.data as abidata
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

def get_projection(fb, l, typat, blist=None, fact=1.0):
    # This function was adapted from the plot_fatbands_typeview function
    e0 = fb.ebands.get_e0('fermie')
    x = np.arange(fb.nkpt)
   
    mybands = range(fb.ebands.mband) if blist is None else blist

    eigens = np.zeros((fb.nkpt, len(mybands)))
    width = np.zeros_like(eigens)

    for symbol in typat:
        for spin in range(fb.nsppol):
            wl_sbk = fb.get_wl_symbol(symbol) * (fact / 2)

            for ib, band in enumerate(mybands):
        
                eigens[:, ib] = fb.ebands.eigens[spin, :, band] - e0
                w = wl_sbk[l, spin, band]
                width[:, ib] += w

    return eigens, width

def plot_relative_projection(eigen, proj1, proj2, ax, x):

    # comptute relative projection and transform to [0,1] interval
    diff = proj1 - proj2
    diff -= np.amin(diff)
    maxproj = np.amax(diff)
    diff = diff/maxproj

    for ib in range(np.shape(diff)[-1]):
        ax0 =ax.scatter(range(x), eigen[:, ib], c=diff[:, ib], marker='o', s=80, cmap='seismic', edgecolor='k', zorder=2, vmin=0, vmax=1)

    return ax0


flist = ["tz2_2o_DS13_FATBANDS.nc", "tz2_2o_DS23_FATBANDS.nc"]

title = ['0GPa', '5GPa']
bands = [np.arange(36,40), np.arange(36,40)]

fig, ax = plt.subplots(1, len(flist), figsize=(5*len(flist), 6))

for j in range(len(flist)):
    # plot electronic bands
    data = abiopen(flist[j])
    data.ebands.set_fermie_to_vbm()
    data.ebands.plot(with_gaps=False, ylims=(-0.4, 0.8), ax=ax[j], show=False, fontsize=16, zorder=1)


    # compute Bi-p projections
    bi_eigen, bi_proj = get_projection(data, l=1, typat=['Bi'], blist=bands[j])
    # compute the sum of Te-p and I-p projections
    tei_eigen, tei_proj = get_projection(data, l=1, typat=['Te', 'I'], blist=bands[j])
    # compute and plot relative projection
    relplot = plot_relative_projection(bi_eigen, bi_proj, tei_proj, x=data.nkpt, ax=ax[j])

    ax[j].set_title(title[j])

    ax[j].set_xticks([0, 17, data.nkpt-1])
    ax[j].set_xticklabels(['0.5H', 'A', '0.5L'], fontsize=12)

divider = make_axes_locatable(ax[1])
cax = divider.append_axes("right", size="7%", pad=0.15)
cbar = plt.colorbar(relplot, cax=cax, ticks=[0, 1])
cbar.ax.set_yticklabels(['(Te+I)$_p$', 'Bi$_p$'])

plt.savefig('tz2_2_relative_projection.png')
plt.show()

