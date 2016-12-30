
import netCDF4 as nc
import matplotlib.pyplot as plt

fname = 'Out/4-2_EP.nc'

with nc.Dataset(fname, 'r') as ds:

    # nomega
    omega = ds.variables['omegase'][...]

    # ntemp
    temperatures = ds.variables['temperatures'][...]

    # nspin, nkpt, nband, nomega, ntemp
    spectral_function = ds.variables['spectral_function_temperature_dependent'][...]


for i, T in enumerate(temperatures):

    # Select one state
    sf = spectral_function[0,0,3,:,i]
    
    plt.plot(omega, sf, label='T={}K'.format(int(T)))
    
plt.legend(loc='upper left')
plt.ylim(0,100)
plt.show()


