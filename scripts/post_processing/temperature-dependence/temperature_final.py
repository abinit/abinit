#!/usr/bin/env python
# Author: Samuel Ponc\'e + Yannick Gillet
# Date: 30/04/2013 -- 11/09/2014 -- 07/08/2015 -- 21/12/2020
# Version: 1.5
# This is the executable script 

#from __future__ import division, print_function
import sys
import os
import copy

try:
  from rf_final import system
  from rf_final import zpm
except ImportError:
  import warnings
  warnings.warn("The rf_final file is missing!")
  raise
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

#############
# Constants #
#############
tol6 = 1E-6
tol8 = 1E-8
Ha2eV = 27.21138386
kb_HaK = 3.1668154267112283e-06

def main():
  start = datetime.now()
  print('Start on %s/%s/%s at %sh%s ' %(start.day,start.month,start.year,start.hour,start.minute))
  
  ######################################################################################
  
  # Interaction with the user
  header = """
    ____  ____       _                                      _                   
   |  _ \|  _ \     | |_ ___ _ __ ___  _ __   ___ _ __ __ _| |_ _   _ _ __ ___  
   | |_) | |_) |____| __/ _ \ '_ ` _ \| '_ \ / _ \ '__/ _` | __| | | | '__/ _ \ 
   |  __/|  __/_____| ||  __/ | | | | | |_) |  __/ | | (_| | |_| |_| | | |  __/ 
   |_|   |_|         \__\___|_| |_| |_| .__/ \___|_|  \__,_|\__|\__,_|_|  \___| 
                                      |_|                              Version 1.5         
  """

  header +=  '\nThis script compute the static/dynamic zero-point motion \n\
  and the temperature dependence of eigenenergies due to electron-phonon interaction.\n\
  The electronic lifetime can also be computed. \n\n\
  WARNING: The first Q-point MUST be the Gamma point.\n'
   
  print(header)

  # Enter the number of cpu on which you want to multi-thread
  user_input = input('Enter the number of cpu on which you want to multi-thread\n')
  nb_cpus = user_input
  try:
    nb_cpus = int(user_input)
  except ValueError:
    raise Exception('The value you enter is not an integer!')
  
  # Type of calculation the user want to perform
  user_input = input('Define the type of calculation you want to perform. Type:\n\
                         1 if you want to run a non-adiabatic AHC calculation\n \
                        2 if you want to run a static AHC calculation\n \
                        3 if you want to run a static AHC calculation without control on active space (not recommended !)\n \
  Note that for 1 & 2 you need _EIGR2D.nc and _GKK.nc files obtained through ABINIT option "ieig2rf 5"\n')
  type = int(user_input)
  
  # Define the output file name
  user_input = input('Enter name of the output file\n')
  output = user_input.strip()
  
  if type <= 0 or type > 3:
    print("Type of the calculation unrecognized !")
    return 1
  
  # Enter the value of the smearing parameter for dynamic AHC
  if type == 1 or type == 2:
    user_input = input('Enter value of the smearing parameter for AHC (in eV)\n')
    smearing = float(user_input)
    smearing = smearing/Ha2eV
  elif type == 3:
    print("The smearing parameter will be the one specified in the abinit input variable elph_imagden\n")
    smearing = None
  else:
    smearing = None
  
  # Enter the value of the smearing parameter for Gaussian broadening
  if type == 1 or type == 2:
    user_input = input('Enter value of the Gaussian broadening for the Eliashberg function and PDOS (in eV)\n')
    gaussian_smearing = float(user_input)
    gaussian_smearing = gaussian_smearing/Ha2eV
  elif type == 3:
    user_input = input('Enter value of the Gaussian broadening for the PDOS (in eV)\n')
    gaussian_smearing = float(user_input)
    gaussian_smearing = gaussian_smearing/Ha2eV
  else:
    gaussian_smearing = None
  
  # Enter the value of the emergy range for the PDOS and Eliashberg calculations
  if type == 1 or type == 2:
    user_input = input('Enter the energy range for the PDOS and Eliashberg calculations (in eV): [e.g. 0 0.5] \n')
    energy_range = user_input.split()
    energy = N.linspace(float(energy_range[0])/Ha2eV,float(energy_range[1])/Ha2eV,500)  
  elif type == 3:
    user_input = input('Enter the energy range for the PDOS (in eV): [e.g. 0 0.5] \n')
    energy_range = user_input.split()
    energy = N.linspace(float(energy_range[0])/Ha2eV,float(energy_range[1])/Ha2eV,500)  
  else:
    energy = None
  
  # Temperature dependence analysis
  user_input = input('Introduce the min temperature, the max temperature and the temperature steps. e.g. 0 200 50 for (0,50,100,150)\n')
  temp_info = user_input.split()
  if(len(temp_info) == 1):
    # If the user specifies only one value, it is the temperature !
    all_temp = [float(temp_info[0])]
  else:
    all_temp = [T for T in N.arange(float(temp_info[0]),float(temp_info[1]),float(temp_info[2]))]
  
  # Lifetime is always activated now !
  lifetime = True
  # Broadening lifetime of the electron
  #user_input = input('Do you want to compute the lifetime of the electrons? [y/n]\n')
  #tmp =user_input.split()[0]
  #if tmp == 'y':
  #  lifetime = True
  #else:
  #  lifetime = False
  
  # Get the nb of Q-points from user 
  user_input = input('Enter the number of Q-points you have\n')
  try:
    nbQ = int(user_input)
  except ValueError:
    raise Exception('The value you enter is not an integer!')

  # SP: update 07/24
  # Get the weight of q-points from the user (not available anymore in _EIGR2D.nc
  user_input = input('Enter the q-point weight for each q-points (in the same order)\n')
  tmp_info = user_input.split()
  wtq = N.zeros((nbQ))
  for ii in N.arange(nbQ):
    wtq[ii] = float(tmp_info[ii])

  # Get the path of the DDB files from user
  DDB_files = []
  for ii in N.arange(nbQ):
    user_input = input('Enter the name of the %s DDB file\n' %ii)
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
    user_input = input('Enter the name of the %s eigq file\n' %ii)
    if len(user_input.split()) != 1:
      raise Exception("You should provide only 1 file")
    else:
      eigq_files.append(user_input.strip(' \t\n\r'))
  
  # Get the path of the EIGR2D files from user
  EIGR2D_files = []
  for ii in N.arange(nbQ):
    user_input = input('Enter the name of the %s EIGR2D file\n' %ii)
    if len(user_input.split()) != 1:
      raise Exception("You should provide only 1 file")
    else:
      EIGR2D_files.append(user_input.strip(' \t\n\r'))
  
  # Get the path of the GKK files from user if dynamical calculation
  if type == 1 or type == 2:
    GKK_files = []
    for ii in N.arange(nbQ):
      user_input = input('Enter the name of the %s GKK file\n' %ii)
      if len(user_input.split()) != 1:
        raise Exception("You should provide only 1 file")
      else:
        GKK_files.append(user_input.strip(' \t\n\r'))
  else:
    GKK_files = []
  
  # Take the EIG at Gamma
  user_input = input('Enter the name of the unperturbed EIG.nc file at Gamma\n')
  if len(user_input.split()) != 1:
    raise Exception("You sould only provide 1 file")
  else:
    eig0 = system(directory='.',filename=user_input.strip(' \t\n\r'))
  
  
  # Read the EIGR2D file at Gamma and save it in ddw_save
  EIGR2D = system(directory='.',filename=EIGR2D_files[0])
  ddw_save = zeros((EIGR2D.nkpt,EIGR2D.nband,3,EIGR2D.natom,3,EIGR2D.natom),dtype=complex)
  ddw_save = copy.deepcopy(EIGR2D.EIG2D_bis)
  
  if type == 1 or type == 2 :
    GKK = system(directory='.',filename=GKK_files[0])
    #print("nkpt = ",GKK.nkpt,"nband = ",GKK.nband,"natom = ",GKK.natom)
    ddw_save2 = zeros((GKK.nkpt,GKK.nband,3,GKK.natom,GKK.nband),dtype=complex)
    ddw_save2 = copy.deepcopy(GKK.GKK_bis)
  else:
    ddw_save2 = None
  
  # Find the degenerate eigenstates
  degen =  zeros((EIGR2D.nsppol,EIGR2D.nkpt,EIGR2D.nband),dtype=int)
  for ispin in N.arange(EIGR2D.nsppol):
    for ikpt in N.arange(EIGR2D.nkpt):
      count = 0
      for iband in N.arange(EIGR2D.nband):
        if iband != EIGR2D.nband-1:
          if N.allclose(eig0.EIG[ispin,ikpt,iband+1], eig0.EIG[ispin,ikpt,iband]):
            degen[ispin,ikpt,iband] = count
          else:
            degen[ispin,ikpt,iband] = count
            count += 1
            continue
        else:
          if N.allclose(eig0.EIG[ispin,ikpt,iband-1], eig0.EIG[ispin,ikpt,iband]):
            degen[ispin,ikpt,iband] = count   
        if iband != 0:
          if N.allclose(eig0.EIG[ispin,ikpt,iband-1], eig0.EIG[ispin,ikpt,iband]):
            degen[ispin,ikpt,iband] = count
        else:
          if N.allclose(eig0.EIG[ispin,ikpt,iband+1], eig0.EIG[ispin,ikpt,iband]):
            degen[ispin,ikpt,iband] = count
  
  # Create the random Q-integration (wtq=1/nqpt):
  #if (abs(EIGR2D.wtq) < tol6):
  #  wtq = N.ones((nbQ))
  #  wtq = wtq*(1.0/nbQ)
  #else:
  #  wtq = N.zeros((nbQ))
  
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
  
  if type == 1 or type == 2:
    arguments = list(zip(nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files,GKK_files))   
    total = zpm(arguments,ddw_save,ddw_save2,nb_cpus,type,\
               all_temp,smearing,eig0_pass,degen,energy,gaussian_smearing)
  elif type == 3:
    arguments = list(zip(nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files))  
    total = zpm(arguments,ddw_save,ddw_save2,nb_cpus,type,\
               all_temp,smearing,eig0_pass,degen,energy,gaussian_smearing)
  total_corr = total.total_corr

  #if (abs(EIGR2D.wtq) > tol6):
  total_wtq = total.total_wtq
  print("Total weigth is ",total_wtq)
  if (total_wtq < 0.9 or total_wtq > 1.1):
    raise Exception("The total weigth is not equal to 1.0. Check that you provide all the q-points.")
  
  
  # Report wall time (before writing final result to be able to include it)
  end = datetime.now()
  print('End on %s/%s/%s at %s h %s ' %(end.day,end.month,end.year,end.hour,end.minute))
  
  runtime = end - start
  print("Runtime: %s seconds (or %s minutes)" %(runtime.seconds,float(runtime.seconds)/60.0))
  
  # Write on a NC files with etsf-io name convention
  ncfile = nc.Dataset(str(output)+'_EP.nc','w')
  # Read dim from first EIGR2D file
  root = nc.Dataset(EIGR2D_files[0],'r')
  # Determine nsppol from reading occ
  nsppol = len(root.variables['occupations'][:,0,0])
  mband = len(root.dimensions['maximum_number_of_bands'])/nsppol
  # Create dimension
  ncfile.createDimension('number_of_atoms',len(root.dimensions['number_of_atoms']))
  ncfile.createDimension('number_of_kpoints',len(root.dimensions['number_of_kpoints']))
  ncfile.createDimension('product_mband_nsppol',len(root.dimensions['maximum_number_of_bands']))
  ncfile.createDimension('cartesian',3)
  ncfile.createDimension('cplex',2)
  ncfile.createDimension('number_of_qpoints',nbQ)
  ncfile.createDimension('number_of_spins',len(root.dimensions['number_of_spins']))
  ncfile.createDimension('max_number_of_states',mband)
  ncfile.createDimension('number_of_temperature',len(all_temp))
  # Create variable
  data = ncfile.createVariable('reduced_coordinates_of_kpoints','d',('number_of_kpoints','cartesian'))
  data[:,:] = root.variables['reduced_coordinates_of_kpoints'][:,:]
  data = ncfile.createVariable('eigenvalues','d',('number_of_spins','number_of_kpoints','max_number_of_states'))
  data[:,:,:] = eig0.EIG
  data = ncfile.createVariable('occupations','i',('number_of_spins','number_of_kpoints','max_number_of_states'))
  data[:,:,:] = root.variables['occupations'][:,:,:]
  data = ncfile.createVariable('primitive_vectors','d',('cartesian','cartesian'))
  data[:,:] = GKK.rprimd
  data = ncfile.createVariable('temperature','d',('number_of_temperature'))
  data[:] = all_temp
  data = ncfile.createVariable('zero_point_motion','d',('number_of_temperature','number_of_spins','number_of_kpoints',\
                              'max_number_of_states','cplex'))
  data[:,:,:,:,0] = total_corr[0,:,:,:,:].real
  data[:,:,:,:,1] = total_corr[0,:,:,:,:].imag
  data = ncfile.createVariable('FAN','d',('number_of_temperature','number_of_spins','number_of_kpoints',\
                              'max_number_of_states','cplex'))
  data[:,:,:,:,0] = total_corr[1,:,:,:,:].real
  data[:,:,:,:,1] = total_corr[1,:,:,:,:].imag
  data = ncfile.createVariable('DW_RIA','d',('number_of_temperature','number_of_spins','number_of_kpoints',\
                              'max_number_of_states','cplex'))
  data[:,:,:,:,0] = total_corr[2,:,:,:,:].real
  data[:,:,:,:,1] = total_corr[2,:,:,:,:].imag
  data = ncfile.createVariable('lifetime','d',('number_of_temperature','number_of_spins','number_of_kpoints',\
                              'max_number_of_states','cplex'))
  data[:,:,:,:,0] = total_corr[3,:,:,:,:].real
  data[:,:,:,:,1] = total_corr[3,:,:,:,:].imag
  # Close the file
  ncfile.close()
  
  # Write Eliashberg file (Does not depend on temperature)
  # Write on a NC files with etsf-io name convention
  
  PHDOS = total_corr[4:len(energy)+4,0,0,0,0]
  if type == 1 or type == 2:
    a2F = total_corr[len(energy)+4:2*len(energy)+4,0,:,:,:]
   
  # Write Eliashberg file (Does not depend on temperature)
  # Write on a NC files with etsf-io name convention
  if type == 1 or type == 2:
    ncfile = nc.Dataset(str(output)+'_g2F.nc','w')
    # Read dim from first EIGR2D file
    root = nc.Dataset(EIGR2D_files[0],'r')
    # Determine nsppol from reading occ
    nsppol = len(root.variables['occupations'][:,0,0])
    mband = len(root.dimensions['maximum_number_of_bands'])/nsppol
    # Create dimension
    ncfile.createDimension('number_of_atoms',len(root.dimensions['number_of_atoms']))
    ncfile.createDimension('number_of_kpoints',len(root.dimensions['number_of_kpoints']))
    ncfile.createDimension('product_mband_nsppol',len(root.dimensions['maximum_number_of_bands']))
    ncfile.createDimension('cartesian',3)
    ncfile.createDimension('cplex',2)
    ncfile.createDimension('number_of_qpoints',nbQ)
    ncfile.createDimension('number_of_spins',len(root.dimensions['number_of_spins']))
    ncfile.createDimension('max_number_of_states',mband)
    ncfile.createDimension('max_number_of_freq',len(energy))
    # Create variable
    data = ncfile.createVariable('reduced_coordinates_of_kpoints','d',('number_of_kpoints','cartesian'))
    data[:,:] = root.variables['reduced_coordinates_of_kpoints'][:,:]
    data = ncfile.createVariable('eigenvalues','d',('number_of_spins','number_of_kpoints','max_number_of_states'))
    data[:,:,:] = eig0.EIG
    data = ncfile.createVariable('occupations','i',('number_of_spins','number_of_kpoints','max_number_of_states'))
    data[:,:,:] = root.variables['occupations'][:,:,:]
    data = ncfile.createVariable('primitive_vectors','d',('cartesian','cartesian'))
    data[:,:] = GKK.rprimd
    data = ncfile.createVariable('Phonon_energy','d',('max_number_of_freq'))
    data[:] = energy
    data = ncfile.createVariable('g2F','d',('number_of_spins','max_number_of_freq','number_of_kpoints',\
                                'max_number_of_states'))
    data[:,:,:,:] = a2F[:,:,:,:].real
    # Close the file
    ncfile.close()
  
  
  # Write PDOS file (Does not depend on temperature)
  # Write on a NC files with etsf-io name convention
  ncfile = nc.Dataset(str(output)+'_PDOS.nc','w')
  # Read dim from first EIGR2D file
  root = nc.Dataset(EIGR2D_files[0],'r')
  # Determine nsppol from reading occ
  nsppol = len(root.variables['occupations'][:,0,0])
  mband = len(root.dimensions['maximum_number_of_bands'])/nsppol
  # Create dimension
  ncfile.createDimension('number_of_atoms',len(root.dimensions['number_of_atoms']))
  ncfile.createDimension('number_of_kpoints',len(root.dimensions['number_of_kpoints']))
  ncfile.createDimension('product_mband_nsppol',len(root.dimensions['maximum_number_of_bands']))
  ncfile.createDimension('cartesian',3)
  ncfile.createDimension('cplex',2)
  ncfile.createDimension('number_of_qpoints',nbQ)
  ncfile.createDimension('number_of_spins',len(root.dimensions['number_of_spins']))
  ncfile.createDimension('max_number_of_states',mband)
  ncfile.createDimension('max_number_of_freq',len(energy))
  # Create variable
  data = ncfile.createVariable('reduced_coordinates_of_kpoints','d',('number_of_kpoints','cartesian'))
  data[:,:] = root.variables['reduced_coordinates_of_kpoints'][:,:]
  data = ncfile.createVariable('eigenvalues','d',('number_of_spins','number_of_kpoints','max_number_of_states'))
  data[:,:,:] = eig0.EIG
  data = ncfile.createVariable('occupations','d',('number_of_spins','number_of_kpoints','max_number_of_states'))
  data[:,:,:] = root.variables['occupations'][:,:,:]
  data = ncfile.createVariable('primitive_vectors','d',('cartesian','cartesian'))
  data[:,:] = GKK.rprimd
  data = ncfile.createVariable('Phonon_energy','d',('max_number_of_freq'))
  data[:] = energy
  data = ncfile.createVariable('phonon_density_of_states','d',('max_number_of_freq'))
  data[:] = PHDOS[:].real
  # Close the file
  ncfile.close()
    
  
  
  # Write the results into the output file
  with open(str(output)+".txt","w") as O:
    O.write("Total correction of the ZPM (eV) for "+str(nbQ)+" Q points\n")
    for ispin in N.arange(nsppol):
      if nsppol > 1:
        if ispin==0: O.write('Spin Up\n')
        else: O.write('Spin down\n')
      for ikpt in N.arange(vkpt):
        O.write('Kpt: '+str(tkpt[ikpt,:])+"\n")
        j = 1
        for ii in (total_corr[0,0,ispin,ikpt,:].real*Ha2eV):
  #       Create a new line every 6 values
          if (j%6 == 0 and j !=0):
            O.write(str(ii)+'\n')
            j += 1
          elif j == vband:
            O.write(str(ii)+'\n')
          else:
            O.write(str(ii)+' ')
            j += 1
    if len(all_temp) > 1:
      O.write("Temperature dependence at Gamma\n")
      for ispin in N.arange(nsppol):
        if nsppol > 1:
          if ispin==0: O.write('Spin Up\n')
          else: O.write('Spin down')
        for iband in N.arange(vband):     
          O.write('Band: '+str(iband)+"\n")
          tt = 0
          for T in all_temp:
            O.write(str(T)+" "+str(total_corr[0,tt,ispin,0,iband].real*Ha2eV)+"\n") 
            tt += 1
    O.write("Fan/DDW contribution at Gamma:\n")
    for ispin in N.arange(nsppol):
      if nsppol > 1:
        if(ispin == 0): O.write('Spin Up\n')
        else: O.write('Spin down\n')
      for iband in N.arange(vband):
        O.write('Band: '+str(iband)+"  FAN: "+str(total_corr[1,0,ispin,0,iband].real*Ha2eV)+"\n")    
        O.write('       '+          "  DDW: "+str(-total_corr[2,0,ispin,0,iband].real*Ha2eV)+"\n")    
        O.write('       '+          "  TOTAL: "+str(total_corr[0,0,ispin,0,iband].real*Ha2eV)+"\n")
    O.write("Runtime: "+str(runtime.seconds)+' seconds (or '+str(float(runtime.seconds)/60.0)+' minutes)')   
  
  # Text file with the PDOS
  with open(str(output)+"_PDOS.txt","w") as O:
    O.write("# Phonon density of state (PDOS) (eV) for "+str(nbQ)+" Q points\n")
    O.write("# Phonon Energy     PDOS\n")
    ii = 0
    for ifreq in energy:
      O.write(str(ifreq*Ha2eV)+' '+str(PHDOS[ii].real)+'\n')  
      ii += 1
    O.write("# Runtime: "+str(runtime.seconds)+' seconds (or '+str(float(runtime.seconds)/60.0)+' minutes)')

  return 0

if __name__ == "__main__":
  sys.exit(main())
