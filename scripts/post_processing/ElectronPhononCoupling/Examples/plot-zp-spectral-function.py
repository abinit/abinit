
import netCDF4 as nc
import matplotlib.pyplot as plt

fname = 'Out/4-1_EP.nc'

with nc.Dataset(fname, 'r') as ds:

    # nomega
    omega = ds.variables['omegase'][...]

    # nspin, nkpt, nband, nomega
    spectral_function = ds.variables['spectral_function'][...]

# Select one state
sf = spectral_function[0,0,3,:]

plt.plot(omega, sf)

plt.show()


