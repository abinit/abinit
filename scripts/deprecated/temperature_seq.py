#!/usr/bin/env python
# Author: Samuel Ponc\'e
# Date: 30/04/2013 -- 11/09/2014
# Version: 1.3
# Script to compute the ZPR

import sys
import os
import copy
try:
  from rf_mods_seq import system
except ImportError:
  import warnings
  warnings.warn("The system module is missing!")
  raise
from rf_mods_seq import zpm
import multiprocessing
from datetime import datetime
try:
  import numpy as N
except ImportError:
  import warnings
  warnings.warn("The numpy module is missing!")
  raise
from numpy import zeros
try:
  import netCDF4 as nc
except ImportError:
  import warnings
  warnings.warn("The netCDF4 module is missing!")
  raise

start = datetime.now()
print 'Start on %s/%s/%s at %sh%s ' %(start.day,start.month,start.year,start.hour,start.minute)

#############
# Constants #
#############
tol6 = 1E-6
tol8 = 1E-8
Ha2eV = 27.21138386
kb_HaK = 3.1668154267112283e-06

######################################################################################

# Interaction with the user
print """
  ____  ____       _                                      _                   
 |  _ \|  _ \     | |_ ___ _ __ ___  _ __   ___ _ __ __ _| |_ _   _ _ __ ___  
 | |_) | |_) |____| __/ _ \ '_ ` _ \| '_ \ / _ \ '__/ _` | __| | | | '__/ _ \ 
 |  __/|  __/_____| ||  __/ | | | | | |_) |  __/ | | (_| | |_| |_| | | |  __/ 
 |_|   |_|         \__\___|_| |_| |_| .__/ \___|_|  \__,_|\__|\__,_|_|  \___| 
                                    |_|                              Version 1.3         
"""
print '\nThis script compute the static/dynamic zero-point motion \n\
and the temperature dependance of eigenenergies due to electron-phonon interaction.\n\
The static electronic lifetime can also be computed. \n\n\
WARNING: The first Q-point MUST be the Gamma point.\n'

# Enter the number of cpu on which you want to multi-thread
user_input = raw_input('Enter the number of cpu on which you want to multi-thread\n')
nb_cpus = user_input
try:
  nb_cpus = int(user_input)
except ValueError:
  raise Exception('The value you enter is not an integer!')

# Type of calculation the user want to perform
user_input = raw_input('Define the type of calculation you want to perform. Type:\n\
                      1 if you want to run a static AHC calculation\n \
                      2 if you want to run a dynamic AHC calculation\n \
                      3 if you want to run a static AHC calculation with control over active space\n\
Note that option 2 and 3 requires _FAN.nc files obtained through ABINIT option "ieig2rf 4\n')
type = N.int(user_input)

# Define the output file name
user_input = raw_input('Enter name of the output file\n')
output = user_input.strip()

# Enter the value of the smearing parameter for dynamic AHC
if (type == 2 or type == 3 ):
  user_input = raw_input('Enter value of the smearing parameter (in eV)\n')
  smearing = N.float(user_input)
  smearing = smearing/Ha2eV
else:
  smearing = None

# Temperature dependence analysis?
user_input = raw_input('Do you want to compute the change of eigenergies with temperature? [y/n]\n')
temperature =user_input.split()[0]
if temperature == 'y':
  temperature = True
  user_input = raw_input('Introduce the starting temperature, max temperature and steps. e.g. 0 2000 100\n')
  temp_info = user_input.split()
else:
  temperature = False
  temp_info = None

# Broadening lifetime of the electron
user_input = raw_input('Do you want to compute the lifetime of the electrons? [y/n]\n')
tmp =user_input.split()[0]
if tmp == 'y':
  lifetime = True
else:
  lifetime = False

# Get the nb of random Q-points from user 
user_input = raw_input('Enter the number of random Q-points you have\n')
try:
  nbQ = int(user_input)
except ValueError:
  raise Exception('The value you enter is not an integer!')

# Get the path of the DDB files from user
DDB_files = []
for ii in N.arange(nbQ):
  user_input = raw_input('Enter the name of the %s DDB file\n' %ii)
  if len(user_input.split()) != 1:
    raise Exception("You should provide only 1 file")
  else: # Append and TRIM the input string with STRIP
    DDB_files.append(user_input.strip(' \t\n\r'))

# Test if the first file is at the Gamma point
DDBtmp = system(directory='.',filename=DDB_files[0])
if N.allclose(DDBtmp.iqpt,[0.0,0.0,0.0]) == False:
  raise Exception('The first Q-point is not Gamma!')

# Get the path of the eigq files from user
eigq_files = []
for ii in N.arange(nbQ):
  user_input = raw_input('Enter the name of the %s eigq file\n' %ii)
  if len(user_input.split()) != 1:
    raise Exception("You should provide only 1 file")
  else:
    eigq_files.append(user_input.strip(' \t\n\r'))

# Get the path of the EIGR2D files from user
EIGR2D_files = []
for ii in N.arange(nbQ):
  user_input = raw_input('Enter the name of the %s EIGR2D file\n' %ii)
  if len(user_input.split()) != 1:
    raise Exception("You should provide only 1 file")
  else:
    EIGR2D_files.append(user_input.strip(' \t\n\r'))

# Get the path of the EIGI2D files from user
if lifetime:
  EIGI2D_files = []
  for ii in N.arange(nbQ):
    user_input = raw_input('Enter the name of the %s EIGI2D file\n' %ii)
    if len(user_input.split()) != 1:
      raise Exception("You should provide only 1 file")
    else:
      EIGI2D_files.append(user_input.strip(' \t\n\r'))

# Get the path of the FAN files from user if dynamical calculation
if (type == 2 or type == 3):
  FAN_files = []
  for ii in N.arange(nbQ):
    user_input = raw_input('Enter the name of the %s FAN file\n' %ii)
    if len(user_input.split()) != 1:
      raise Exception("You should provide only 1 file")
    else:
      FAN_files.append(user_input.strip(' \t\n\r'))

# Take the EIG at Gamma
user_input = raw_input('Enter the name of the unperturbed EIG.nc file at Gamma\n')
if len(user_input.split()) != 1:
  raise Exception("You sould only provide 1 file")
else:
  eig0 = system(directory='.',filename=user_input.strip(' \t\n\r'))


# Read the EIGR2D file at Gamma and save it in ddw_save
EIGR2D = system(directory='.',filename=EIGR2D_files[0])
ddw_save = zeros((EIGR2D.nkpt,EIGR2D.nband,3,EIGR2D.natom,3,EIGR2D.natom),dtype=complex)
ddw_save = copy.deepcopy(EIGR2D.EIG2D)

if (type == 2 or type == 3):
  FAN = system(directory='.',filename=FAN_files[0])
  ddw_save2 = zeros((FAN.nkpt,FAN.nband,3,FAN.natom,3,FAN.natom,FAN.nband),dtype=complex)
  ddw_save2 = copy.deepcopy(FAN.FAN)
if (type == 1):
# We put dummy argument
  ddw_save2 = 0.0

# Find the degenerate eigenstates
degen =  zeros((EIGR2D.nkpt,EIGR2D.nband),dtype=int)
for ikpt in N.arange(EIGR2D.nkpt):
  count = 0
  for iband in N.arange(EIGR2D.nband):
    if iband != EIGR2D.nband-1:
      if N.allclose(eig0.EIG[0,ikpt,iband+1], eig0.EIG[0,ikpt,iband]):
        degen[ikpt,iband] = count
      else:
        degen[ikpt,iband] = count
        count += 1
        continue
    else:
      if N.allclose(eig0.EIG[0,ikpt,iband-1], eig0.EIG[0,ikpt,iband]):
        degen[ikpt,iband] = count   
    if iband != 0:
      if N.allclose(eig0.EIG[0,ikpt,iband-1], eig0.EIG[0,ikpt,iband]):
        degen[ikpt,iband] = count
    else:
      if N.allclose(eig0.EIG[0,ikpt,iband+1], eig0.EIG[0,ikpt,iband]):
        degen[ikpt,iband] = count

# Create the random Q-integration (wtq=1/nqpt):
if (EIGR2D.wtq == 0):
  wtq = N.ones((nbQ))
  wtq = wtq*(1.0/nbQ)
else:
  wtq = N.zeros((nbQ))

#DBSP
#wtq = N.ones((nbQ))
#END
nbqpt = N.arange(nbQ)
  
# Compute phonon freq. and eigenvector for each Q-point 
# from each DDB (1 qpt per DDB file)
vkpt = EIGR2D.nkpt
vband = EIGR2D.nband
tkpt = zeros((EIGR2D.nkpt,3))
tkpt = EIGR2D.kpt[:,:]
eig0_pass = copy.deepcopy(eig0.EIG)

if temperature:
  temp_info = N.arange(N.float(temp_info[0]),N.float(temp_info[1]),N.float(temp_info[2]))

if (type == 1):
  total = zpm(zip(nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files),ddw_save,ddw_save2,nb_cpus,type,temperature,\
             temp_info,smearing,eig0_pass,degen,lifetime=False)
  if lifetime:
    print "Now compute broadening ..."
    broadening = zpm(zip(nbqpt,wtq,eigq_files,DDB_files,EIGI2D_files),ddw_save,ddw_save2,nb_cpus,type,temperature,\
            temp_info,smearing,eig0_pass,degen,lifetime)
if (type == 2):
  if ((FAN.nkpt*FAN.nband*3*FAN.natom*3*FAN.natom*FAN.nband)*8.0)/(1024**2) < 1000: # If FAN file is smaller than 1000MB
    total = zpm(zip(nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files,FAN_files),ddw_save,ddw_save2,nb_cpus,type,temperature,\
             temp_info,smearing,eig0_pass,degen,lifetime=False)
  else:  # Use sequential version without multiprocess threading ==> otherwise it breaks down
    total = zpm(zip(nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files,FAN_files),ddw_save,ddw_save2,1,type,temperature,\
             temp_info,smearing,eig0_pass,degen,lifetime=False)
  if lifetime:
    #raise Exception("Dynamical lifetime is not yet implemented.")
    print "WARNING:"
    print "Dynamical lifetime is not yet implemented...proceed with static lifetime"
    print "Now compute broadening ..."
    broadening = zpm(zip(nbqpt,wtq,eigq_files,DDB_files,EIGI2D_files),ddw_save,ddw_save2,nb_cpus,type,temperature,\
            temp_info,smearing,eig0_pass,degen,lifetime)

if (type == 3):
  if ((FAN.nkpt*FAN.nband*3*FAN.natom*3*FAN.natom*FAN.nband)*8.0)/(1024**2) < 1000: # If FAN file is smaller than 1000MB
    total = zpm(zip(nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files,FAN_files),ddw_save,ddw_save2,nb_cpus,type,temperature,\
             temp_info,smearing,eig0_pass,degen,lifetime=False)
  else: # Use sequential version without multiprocess threading ==> otherwise it breaks down
    total = zpm(zip(nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files,FAN_files),ddw_save,ddw_save2,1,type,temperature,\
             temp_info,smearing,eig0_pass,degen,lifetime=False)
  if lifetime:
    print "Now compute broadening ..."
    broadening = zpm(zip(nbqpt,wtq,eigq_files,DDB_files,EIGI2D_files),ddw_save,ddw_save2,nb_cpus,type,temperature,\
            temp_info,smearing,eig0_pass,degen,lifetime)

total_corr = total.total_corr 

if lifetime:
  brd_total = broadening.broadening  

if (EIGR2D.wtq != 0):
  total_wtq = total.total_wtq
  print "Total weigth is ",total_wtq
  if (total_wtq < 0.9 or total_wtq > 1.1):
    raise Exception("The total weigth is not equal to 1.0. Check that you provide all the q-points.")


# Report wall time (before writing final result to be able to include it)
end = datetime.now()
print 'End on %s/%s/%s at %s h %s ' %(end.day,end.month,end.year,end.hour,end.minute)

runtime = end - start
print "Runtime: %s seconds (or %s minutes)" %(runtime.seconds,float(runtime.seconds)/60.0)

if temperature:
  if lifetime:
    # Write on a NC files with etsf-io name convention
    ncfile = nc.Dataset(str(output)+'_EP.nc','w')
    # Read dim from first EIGR2D file
    root = nc.Dataset(EIGR2D_files[0],'r')
    # Determine nsppol from reading occ
    nsppol = len(root.variables['occupations'][:,0,0])
    if nsppol > 1:
      print "WARNING: nsppol > 1 has not been tested."
    mband = len(root.dimensions['product_mband_nsppol'])/nsppol
    # Create dimension
    ncfile.createDimension('number_of_atoms',len(root.dimensions['number_of_atoms']))
    ncfile.createDimension('number_of_kpoints',len(root.dimensions['number_of_kpoints']))
    ncfile.createDimension('product_mband_nsppol',len(root.dimensions['product_mband_nsppol']))
    ncfile.createDimension('cartesian',3)
    ncfile.createDimension('cplex',2)
    ncfile.createDimension('number_of_qpoints',nbQ)
    ncfile.createDimension('number_of_spins',len(root.dimensions['number_of_spins']))
    ncfile.createDimension('max_number_of_states',mband)
    ncfile.createDimension('number_of_temperature',len(temp_info))
    # Create variable
    data = ncfile.createVariable('reduced_coordinates_of_kpoints','d',('number_of_kpoints','cartesian'))
    data[:,:] = root.variables['reduced_coordinates_of_kpoints'][:,:]
    data = ncfile.createVariable('eigenvalues','d',('number_of_spins','number_of_kpoints','max_number_of_states'))
    data[:,:,:] = root.variables['eigenvalues'][:,:,:]
    data = ncfile.createVariable('occupations','i',('number_of_spins','number_of_kpoints','max_number_of_states'))
    data[:,:,:] = root.variables['occupations'][:,:,:]
    data = ncfile.createVariable('primitive_vectors','d',('cartesian','cartesian'))
    data[:,:] = root.variables['primitive_vectors'][:,:]
    data = ncfile.createVariable('temperature','d',('number_of_temperature'))
    data[:] = temp_info
    data = ncfile.createVariable('zero_point_motion','d',('number_of_spins','number_of_temperature','number_of_kpoints',\
                                'max_number_of_states','cplex'))
    data[0,:,:,:,0] = total_corr[0,:,:,:].real
    data[0,:,:,:,1] = brd_total[:,:,:].real 
    # Close the file
    ncfile.close()
  else:
    # Write on a NC files with etsf-io name convention
    ncfile = nc.Dataset(str(output)+'_EP.nc','w')
    # Read dim from first EIGR2D file
    root = nc.Dataset(EIGR2D_files[0],'r')
    # Determine nsppol from reading occ
    nsppol = len(root.variables['occupations'][:,0,0])
    if nsppol > 1:
      print "WARNING: nsppol > 1 has not been tested."
    mband = len(root.dimensions['product_mband_nsppol'])/nsppol
    # Create dimension
    ncfile.createDimension('number_of_atoms',len(root.dimensions['number_of_atoms']))
    ncfile.createDimension('number_of_kpoints',len(root.dimensions['number_of_kpoints']))
    ncfile.createDimension('product_mband_nsppol',len(root.dimensions['product_mband_nsppol']))
    ncfile.createDimension('cartesian',3)
    ncfile.createDimension('cplex',2)
    ncfile.createDimension('number_of_qpoints',nbQ)
    ncfile.createDimension('number_of_spins',len(root.dimensions['number_of_spins']))
    ncfile.createDimension('max_number_of_states',mband)
    ncfile.createDimension('number_of_temperature',len(temp_info))
    # Create variable
    data = ncfile.createVariable('reduced_coordinates_of_kpoints','d',('number_of_kpoints','cartesian'))
    data[:,:] = root.variables['reduced_coordinates_of_kpoints'][:,:]
    data = ncfile.createVariable('eigenvalues','d',('number_of_spins','number_of_kpoints','max_number_of_states'))
    data[:,:,:] = root.variables['eigenvalues'][:,:,:]
    data = ncfile.createVariable('occupations','i',('number_of_spins','number_of_kpoints','max_number_of_states'))
    data[:,:,:] = root.variables['occupations'][:,:,:]
    data = ncfile.createVariable('primitive_vectors','d',('cartesian','cartesian'))
    data[:,:] = root.variables['primitive_vectors'][:,:]
    data = ncfile.createVariable('temperature','d',('number_of_temperature'))
    data[:] = temp_info
    data = ncfile.createVariable('zero_point_motion','d',('number_of_spins','number_of_temperature','number_of_kpoints',\
                                'max_number_of_states','cplex'))
    data[0,:,:,:,0] = total_corr[0,:,:,:].real
    data[0,:,:,:,1] = 0.0
    # Close the file
    ncfile.close()
else:
  if lifetime:
    # Write on a NC files with etsf-io name convention
    ncfile = nc.Dataset(str(output)+'_EP.nc','w')
    # Read dim from first EIGR2D file
    root = nc.Dataset(EIGR2D_files[0],'r')
    # Determine nsppol from reading occ
    nsppol = len(root.variables['occupations'][:,0,0])
    if nsppol > 1:
      print "WARNING: nsppol > 1 has not been tested."
    mband = len(root.dimensions['product_mband_nsppol'])/nsppol
    # Create dimension
    ncfile.createDimension('number_of_atoms',len(root.dimensions['number_of_atoms']))
    ncfile.createDimension('number_of_kpoints',len(root.dimensions['number_of_kpoints']))
    ncfile.createDimension('product_mband_nsppol',len(root.dimensions['product_mband_nsppol']))
    ncfile.createDimension('cartesian',3)
    ncfile.createDimension('cplex',2)
    ncfile.createDimension('number_of_qpoints',nbQ)
    ncfile.createDimension('number_of_spins',len(root.dimensions['number_of_spins']))
    ncfile.createDimension('max_number_of_states',mband)
    ncfile.createDimension('number_of_temperature',1)
    # Create variable
    data = ncfile.createVariable('reduced_coordinates_of_kpoints','d',('number_of_kpoints','cartesian'))
    data[:,:] = root.variables['reduced_coordinates_of_kpoints'][:,:]
    data = ncfile.createVariable('eigenvalues','d',('number_of_spins','number_of_kpoints','max_number_of_states'))
    data[:,:,:] = root.variables['eigenvalues'][:,:,:]
    data = ncfile.createVariable('occupations','i',('number_of_spins','number_of_kpoints','max_number_of_states'))
    data[:,:,:] = root.variables['occupations'][:,:,:]
    data = ncfile.createVariable('primitive_vectors','d',('cartesian','cartesian'))
    data[:,:] = root.variables['primitive_vectors'][:,:]
    data = ncfile.createVariable('temperature','d',('number_of_temperature'))
    data[:] = 0.0
    data = ncfile.createVariable('zero_point_motion','d',('number_of_spins','number_of_temperature','number_of_kpoints',\
                                'max_number_of_states','cplex'))
    data[0,0,:,:,0] = total_corr[0,:,:].real
    data[0,0,:,:,1] = brd_total[:,:].real 
    # Close the file
    ncfile.close() 
  else: 
    # Write on a NC files with etsf-io name convention
    ncfile = nc.Dataset(str(output)+'_EP.nc','w')
    # Read dim from first EIGR2D file
    root = nc.Dataset(EIGR2D_files[0],'r')
    # Determine nsppol from reading occ
    nsppol = len(root.variables['occupations'][:,0,0])
    if nsppol > 1:
      print "WARNING: nsppol > 1 has not been tested."
    mband = len(root.dimensions['product_mband_nsppol'])/nsppol
    # Create dimension
    ncfile.createDimension('number_of_atoms',len(root.dimensions['number_of_atoms']))
    ncfile.createDimension('number_of_kpoints',len(root.dimensions['number_of_kpoints']))
    ncfile.createDimension('product_mband_nsppol',len(root.dimensions['product_mband_nsppol']))
    ncfile.createDimension('cartesian',3)
    ncfile.createDimension('cplex',2)
    ncfile.createDimension('number_of_qpoints',nbQ)
    ncfile.createDimension('number_of_spins',len(root.dimensions['number_of_spins']))
    ncfile.createDimension('max_number_of_states',mband)
    ncfile.createDimension('number_of_temperature',1)
    # Create variable
    data = ncfile.createVariable('reduced_coordinates_of_kpoints','d',('number_of_kpoints','cartesian'))
    data[:,:] = root.variables['reduced_coordinates_of_kpoints'][:,:]
    data = ncfile.createVariable('eigenvalues','d',('number_of_spins','number_of_kpoints','max_number_of_states'))
    data[:,:,:] = root.variables['eigenvalues'][:,:,:]
    data = ncfile.createVariable('occupations','i',('number_of_spins','number_of_kpoints','max_number_of_states'))
    data[:,:,:] = root.variables['occupations'][:,:,:]
    data = ncfile.createVariable('primitive_vectors','d',('cartesian','cartesian'))
    data[:,:] = root.variables['primitive_vectors'][:,:]
    data = ncfile.createVariable('temperature','d',('number_of_temperature'))
    data[:] = 0.0
    data = ncfile.createVariable('zero_point_motion','d',('number_of_spins','number_of_temperature','number_of_kpoints',\
                               'max_number_of_states','cplex'))
    data[0,0,:,:,0] = total_corr[0,:,:].real
    data[0,0,:,:,1] = 0.0
    # Close the file
    ncfile.close()

# Write the results into the output file
if temperature:
  with open(str(output)+".txt","w") as O:
    O.write("Total correction of the ZPM (eV) for "+str(nbQ)+" Q points\n")
    for ikpt in N.arange(vkpt):
      O.write('Kpt: '+str(tkpt[ikpt,:])+"\n")
      j = 1
      for ii in (total_corr[0,0,ikpt,:].real*Ha2eV):
  #     Create a new line every 6 values
        if (j%6 == 0 and j !=0):
          O.write(str(ii)+'\n')
          j += 1
        elif j == vband:
          O.write(str(ii)+'\n')
        else:
          O.write(str(ii)+' ')
          j += 1
    O.write("Temperature dependence at Gamma\n")
    for iband in N.arange(vband):     
      O.write('Band: '+str(iband)+"\n")
      tt = 0
      for T in temp_info:
        O.write(str(T)+" "+str(total_corr[0,tt,0,iband].real*Ha2eV)+"\n") 
	tt += 1
    O.write("Fan/DDW contribution at Gamma:\n")
    for iband in N.arange(vband):
      O.write('Band: '+str(iband)+"  FAN: "+str(total_corr[1,0,0,iband].real*Ha2eV)+"\n")    
      O.write('       '+          "  DDW: "+str(-total_corr[2,0,0,iband].real*Ha2eV)+"\n")    
      O.write('       '+          "  TOTAL: "+str(total_corr[0,0,0,iband].real*Ha2eV)+"\n")
    O.write("Runtime: "+str(runtime.seconds)+' seconds (or '+str(float(runtime.seconds)/60.0)+' minutes)')   
  if lifetime:
    with open(str(output)+"_BRD.txt","w") as O:
      O.write("Total correction of the ZPM (eV) for "+str(nbQ)+" Q points\n")
      for ikpt in N.arange(vkpt):
        O.write('Kpt: '+str(tkpt[ikpt,:])+"\n")
        j = 1
        for ii in (brd_total[0,ikpt,:].real*Ha2eV):
    #     Create a new line every 6 values
          if (j%6 == 0 and j !=0):
            O.write(str(ii)+'\n')
            j += 1
          elif j == vband:
            O.write(str(ii)+'\n')
          else:
            O.write(str(ii)+' ')
            j += 1
      O.write("Temperature dependence at Gamma\n")
      for iband in N.arange(vband):
        O.write('Band: '+str(iband)+"\n")
        tt = 0
        for T in temp_info: 
          O.write(str(T)+" "+str(brd_total[tt,0,iband].real*Ha2eV)+"\n")
          tt += 1
      O.write("Runtime: "+str(runtime.seconds)+' seconds (or '+str(float(runtime.seconds)/60.0)+' minutes)')
else:
  with open(str(output)+".txt","w") as O:
    O.write("Total correction of the ZPM (eV) for "+str(nbQ)+" Q points\n")
    for ikpt in N.arange(vkpt):
      O.write('Kpt: '+str(tkpt[ikpt,:])+"\n")
      j = 1
      for ii in (total_corr[0,ikpt,:].real*Ha2eV):
  #     Create a new line every 6 values
        if (j%6 == 0 and j !=0):
          O.write(str(ii)+'\n')
          j += 1
        elif j == vband:
          O.write(str(ii)+'\n')
        else:
          O.write(str(ii)+' ')
          j += 1
    O.write("Fan/DDW contribution at Gamma:\n")
    for iband in N.arange(vband):
      O.write('Band: '+str(iband)+"  FAN: "+str(total_corr[1,0,iband].real*Ha2eV)+"\n")    
      O.write('       '+           "  DDW: "+str(-total_corr[2,0,iband].real*Ha2eV)+"\n")    
      O.write('       '+           "  TOTAL: "+str(total_corr[0,0,iband].real*Ha2eV)+"\n")    
    O.write("Runtime: "+str(runtime.seconds)+' seconds (or '+str(float(runtime.seconds)/60.0)+' minutes)')   
  if lifetime:
    with open(str(output)+"_BRD.txt","w") as O:
      O.write("Total correction of the ZPM (eV) for "+str(nbQ)+" Q points\n")
      for ikpt in N.arange(vkpt):
        O.write('Kpt: '+str(tkpt[ikpt,:])+"\n")
        j = 1
        for ii in (brd_total[ikpt,:].real*Ha2eV):
#         Create a new line every 6 values
          if (j%6 == 0 and j !=0):
            O.write(str(ii)+'\n')
            j += 1
          elif j == vband:
            O.write(str(ii)+'\n')
          else:
            O.write(str(ii)+' ')
            j += 1
      O.write("Runtime: "+str(runtime.seconds)+' seconds (or '+str(float(runtime.seconds)/60.0)+' minutes)')






