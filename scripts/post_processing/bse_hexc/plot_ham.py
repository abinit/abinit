#!/bin/bash

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import copy
from hexc import Hexc

if __name__ == '__main__':
    ham2 = Hexc.from_file('t31o_DS3_HEXC.nc')
    ham3 = Hexc.from_file('t31o_DS5_HEXC.nc')
    ham4 = Hexc.from_file('t32o_DS1_HEXC_I.nc')
    ham5 = Hexc.from_file('t33o_DS1_HEXC_I.nc')
    ham6 = Hexc.from_file('t34o_DS1_HEXC_I.nc')

    list_ikpt2 = None
    list_ikpt3 = None
    list_ikpt4 = None
    list_ikpt5 = None
    list_ikpt6 = None

    #yval = 0.06
    #zval = 0.065

    #list_ikpt2 = ham2.get_filtered_kpt(yval=yval,zval=zval)
    #list_ikpt3 = ham3.get_filtered_kpt(yval=yval,zval=zval)
    #list_ikpt4 = ham4.get_filtered_kpt(yval=yval,zval=zval)
    #list_ikpt5 = ham5.get_filtered_kpt(yval=yval,zval=zval)
    #list_ikpt6 = ham6.get_filtered_kpt(yval=yval,zval=zval)
    
    ham2.ham[:,:] = ham2.ham[:,:] - np.diag(ham2.diag[:])
    ham3.ham[:,:] = ham3.ham[:,:] - np.diag(ham3.diag[:])

    icval = 4
    ivval = 1
    icpval = 4
    ivpval = 1

    ham2.plot(ic=icval,iv=ivval,icp=icpval,ivp=ivpval,list_ikpt=list_ikpt2,title='Coarse')
    ham3.plot(ic=icval,iv=ivval,icp=icpval,ivp=ivpval,list_ikpt=list_ikpt3,title='Dense')
    ham4.plot(ic=icval,iv=ivval,icp=icpval,ivp=ivpval,list_ikpt=list_ikpt4,title='M1')
    ham5.plot(ic=icval,iv=ivval,icp=icpval,ivp=ivpval,list_ikpt=list_ikpt5,title='M2')
    ham6.plot(ic=icval,iv=ivval,icp=icpval,ivp=ivpval,list_ikpt=list_ikpt6,title='M3')

    plt.show()





