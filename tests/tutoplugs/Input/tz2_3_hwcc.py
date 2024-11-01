#!/usr/bin/env python
import os
import matplotlib.pyplot as plt
import z2pack

result_0 = z2pack.io.load('results_tz2_3/BiTeI_0.msgpack')
result_1 = z2pack.io.load('results_tz2_3/BiTeI_1.msgpack')
result_2 = z2pack.io.load('results_tz2_3/BiTeI_2.msgpack')
result_3 = z2pack.io.load('results_tz2_3/BiTeI_3.msgpack')

# plotting
fig, ax = plt.subplots(2, 2, figsize=(10,10), constrained_layout=True)

# plot styling
fs = 15
for i in range(2):
    for j in range(2):
        ax[i,j].set_xlabel(r'$k_x$', fontsize=fs, labelpad=-5)
        ax[i,j].set_xticks([0, 1])
        ax[i,j].set_xticklabels([r'$0$', r'$\pi/a_x$'])
ax[0,0].set_ylabel(r'$\bar{y}$', rotation='horizontal', fontsize=fs)
ax[0,0].set_title(
    r'$0 GPa, k_z=0.0, \Delta={}$'.format(z2pack.invariant.z2(result_0)),
    fontsize=fs
)
ax[0,1].set_title(
    r'$ 0 GPa, k_z=0.5, \Delta={}$'.format(z2pack.invariant.z2(result_1)),
    fontsize=fs
)
ax[1,0].set_ylabel(r'$\bar{y}$', rotation='horizontal', fontsize=fs)
ax[1,0].set_title(
    r'$5 GPa, k_z=0.0, \Delta={}$'.format(z2pack.invariant.z2(result_2)),
    fontsize=fs
)
ax[1,1].set_title(
    r'$ 5 GPa, k_z=0.5, \Delta={}$'.format(z2pack.invariant.z2(result_3)),
    fontsize=fs
)
# plotting the WCC evolution
z2pack.plot.wcc(result_0, axis=ax[0,0])
z2pack.plot.wcc(result_1, axis=ax[0,1])
z2pack.plot.wcc(result_2, axis=ax[1,0])
z2pack.plot.wcc(result_3, axis=ax[1,1])
#z2pack.plot.chern(result_1, axis=ax[0,2])
#z2pack.plot.chern(result_3, axis=ax[1,2])


#if not os.path.isdir('plots'):
#    os.mkdir('plots')

plt.savefig('tz2_3_hwcc.png', bbox_inches='tight')
plt.show()
