
import netCDF4 as nc
import matplotlib.pyplot as plt

fname = 'Out/4-1_EP.nc'

with nc.Dataset(fname, 'r') as ds:

    # nomega
    omega = ds.variables['omegase'][...]

    # nspin, nkpt, nband, nomega
    self_energy = ds.variables['self_energy'][...]

# Select one state
re = self_energy[0,0,3,:,0]
im = self_energy[0,0,3,:,1]

plt.plot(omega, re)
plt.plot(omega, im)

plt.show()


