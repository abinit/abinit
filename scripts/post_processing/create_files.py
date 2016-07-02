#!/usr/bin/env python
# Author: Samuel Ponc\'e
# Date: 24/04/2013
# Script create the files file to compute the ZPR

import numpy as np

nbQ = 2500
files_name = 'temperature.files'

with open(files_name,'w') as O:
  O.write("1\n")
  O.write("ZPM_para_1\n")
  O.write("0.1\n")
  O.write("n\n")
  O.write(str(nbQ)+"\n")
  for ii in np.arange(nbQ):
    O.write("Q_"+str(ii)+"/abo_DS3_DDB\n")
  for ii in np.arange(nbQ):
    O.write("Q_"+str(ii)+"/abo_DS2_EIG.nc\n")
  for ii in np.arange(nbQ):
    O.write("Q_"+str(ii)+"/abo_DS3_EIGR2D\n")
  O.write("Q_0/abo_DS1_EIG.nc")
