#!/usr/bin/env python
# Author: Samuel Ponc\'e
# Date: 30/04/2013
# Script to compute the spectral function

import sys
import os
from rf_mods import system
import multiprocessing
from datetime import datetime
try:
  import numpy as N
except ImportError:
  import warnings
  warnings.warn("The numpy module is missing!")
  raise
try:
  import netCDF4 as nc
except ImportError:
  import warnings
  warnings.warn("The netCDF4 module is missing!")
  raise
#import matplotlib.pyplot as plt


start = datetime.now()
print 'Start on %s/%s/%s at %sh%s ' %(start.day,start.month,start.year,start.hour,start.minute)

#############
# Constants #
#############
# If you want to hardcode the weight of the k-points you can do it here:
# wtq = [0.00, 0.125, 0.5,0.375]

# The numb of cpu used is hardcoded for now 
nb_cpus = 2

tol6 = 1E-6
tol8 = 1E-8
Ha2eV = 27.21138386
kb_HaK = 3.1668154267112283e-06

##############
# Definition #
##############

def FanDDW(arguments):
  wtq,eigq_files,DDB_files,EIGR2D_files,FAN_files = arguments
  DDB = system()
  FANterm = system()
  EIGR2D = system()
  eigq = system()
  DDB.__init__(directory='.',filename=DDB_files)
  tot_corr = N.zeros((len(freq)),dtype=complex)
  self_energy = N.zeros((len(freq)),dtype=complex)
  spectral = N.zeros((len(freq)))


# Calcul of gprimd from rprimd
  rprimd = DDB.rprim*DDB.acell
  gprimd = N.linalg.inv(N.matrix(rprimd))

# Transform from 2nd-order matrix (non-cartesian coordinates, 
# masses not included, asr not included ) from DDB to
# dynamical matrix, in cartesian coordinates, asr not imposed.
  IFC_cart = N.zeros((3,DDB.natom,3,DDB.natom),dtype=complex)
  for ii in N.arange(DDB.natom):
    for jj in N.arange(DDB.natom):
      for dir1 in N.arange(3):
        for dir2 in N.arange(3):
          for dir3 in N.arange(3):
            for dir4 in N.arange(3):
              IFC_cart[dir1,ii,dir2,jj] += gprimd[dir1,dir3]*DDB.IFC[dir3,ii,dir4,jj] \
            *gprimd[dir2,dir4]

# Reduce the 4 dimensional IFC_cart matrice to 2 dimensional Dynamical matrice.
  ipert1 = 0
  Dyn_mat = N.zeros((3*DDB.natom,3*DDB.natom),dtype=complex)
  while ipert1 < 3*DDB.natom:
    for ii in N.arange(DDB.natom):
      for dir1 in N.arange(3):
        ipert2 = 0
        while ipert2 < 3*DDB.natom:
          for jj in N.arange(DDB.natom):
            for dir2 in N.arange(3):
              Dyn_mat[ipert1,ipert2] = IFC_cart[dir1,ii,dir2,jj]
              ipert2 += 1
        ipert1 += 1

# Hermitianize the dynamical matrix
  dynmat = N.matrix(Dyn_mat)
  dynmat = 0.5*(dynmat + dynmat.transpose().conjugate())

# Solve the eigenvalue problem with linear algebra (Diagonalize the matrix)
  [eigval,eigvect]=N.linalg.eigh(Dyn_mat)

# Orthonormality relation 
  eigvect = (eigvect)*N.sqrt(5.4857990965007152E-4/float(DDB.amu[0]))

# Phonon frequency (5.4857990946E-4 = 1 au of electron mass)
  omega = N.sqrt((eigval*5.4857990965007152E-4)/float(DDB.amu[0]))

# Now read the EIGq, EIGR2D and FAN   
  eigq.__init__(directory='.',filename=eigq_files)
  EIGR2D.__init__(directory='.',filename=EIGR2D_files)
  FANterm.__init__(directory='.',filename=FAN_files)

# Compute the displacement = eigenvectors of the DDB. 
# Due to metric problem in reduce coordinate we have to work in cartesian
# but then go back to reduce because our EIGR2D matrix elements are in reduced coord.
  displ_FAN =  N.zeros((3,3),dtype=complex)
  displ_DDW =  N.zeros((3,3),dtype=complex)
  fan_add = N.zeros((len(freq)),dtype=complex)
  fan_corr = N.zeros((),dtype=complex)
  ddw_corr = N.zeros((),dtype=complex)

  for imode in N.arange(3*EIGR2D.natom): #Loop on perturbation (6 for 2 atoms)
    if omega[imode].real > tol6:
      for iatom1 in N.arange(EIGR2D.natom):
        for iatom2 in N.arange(EIGR2D.natom):
          for idir1 in N.arange(0,3):
            for idir2 in N.arange(0,3):
              displ_FAN[idir1,idir2] = eigvect[3*iatom2+idir2,imode].conj()\
                 *eigvect[3*iatom1+idir1,imode]/(2.0*omega[imode].real)
              displ_DDW[idir1,idir2] = (eigvect[3*iatom2+idir2,imode].conj()\
                 *eigvect[3*iatom2+idir1,imode]+eigvect[3*iatom1+idir2,imode].conj()\
                 *eigvect[3*iatom1+idir1,imode])/(4.0*omega[imode].real)
          # Now switch to reduced coordinates in 2 steps (more efficient)
          tmp_displ_FAN = N.zeros((3,3),dtype=complex)
          tmp_displ_DDW = N.zeros((3,3),dtype=complex)
          for idir1 in N.arange(3):
            for idir2 in N.arange(3):
              tmp_displ_FAN[:,idir1] = tmp_displ_FAN[:,idir1]+displ_FAN[:,idir2]*gprimd[idir2,idir1]
              tmp_displ_DDW[:,idir1] = tmp_displ_DDW[:,idir1]+displ_DDW[:,idir2]*gprimd[idir2,idir1]
          displ_red_FAN = N.zeros((3,3),dtype=complex)
          displ_red_DDW = N.zeros((3,3),dtype=complex)
          for idir1 in N.arange(3):
            for idir2 in N.arange(3):
              displ_red_FAN[idir1,:] = displ_red_FAN[idir1,:] + tmp_displ_FAN[idir2,:]*gprimd[idir2,idir1]
              displ_red_DDW[idir1,:] = displ_red_DDW[idir1,:] + tmp_displ_DDW[idir2,:]*gprimd[idir2,idir1]
          # Now compute the T=0 shift due to this q point
          for idir1 in N.arange(3):
            for idir2 in N.arange(3):
              fan_corr += EIGR2D.EIG2D[kpt,band-1,idir1,iatom1,idir2,iatom2]*\
                  displ_red_FAN[idir1,idir2]
              ddw_corr += ddw_save[idir1,iatom1,idir2,iatom2]*\
                  displ_red_DDW[idir1,idir2]
      if temperature < tol6:
        bose = 0
      else:
        bose = 1.0/(N.exp(omega[imode].real/(kb_HaK*temperature))-1)
      
      fan_corr = fan_corr*(2*bose+1.0)
      ddw_corr = ddw_corr*(2*bose+1.0)
      for jband in N.arange(EIGR2D.nband):
        index = 0
        for ifreq in freq:
          delta_E = ifreq - eigq.EIG[0,ikpt,jband] + smearing*1j
          fan_add[index] = fan_add[index] + FANterm.FAN[kpt,band-1,imode,jband]*(\
                                  (bose+0.5)*(2*delta_E/(delta_E**2-(omega[imode].real)**2)) \
                                - (1-EIGR2D.occ[iband])*(omega[imode].real/(delta_E**2-(omega[imode].real)**2))\
                                 -(bose+0.5)*2/delta_E)/(2.0*omega[imode].real)
          index += 1

  tot_corr[:] = (fan_corr+ddw_corr+fan_add[:])*wtq
  return tot_corr
#####################
# End of definitions
######################################################################################

# Interaction with the user
print '\n##################################################'
print '# Spectral function of the dynamical self-energy #'
print '##################################################'
print '\nThis script compute the zero-point motion and the temperature dependence \n\
of eigenenergies due to electron-phonon interaction. This script can \n\
only compute Q-points with the same weight. If you want symmetry you must hack the script.\n\
WARNING: The first Q-point MUST be the Gamma point\n'

# Define the output file name
user_input = raw_input('Enter name of the output file\n')
output = user_input

# Get the value of the smearing parameter (in eV)
user_input = raw_input('Enter value of the smearing parameter (in eV)\n')
smearing = N.float(user_input)
smearing = smearing/Ha2eV

# Frequency range?
user_input = raw_input('What is the range of frequencies you want to compute your spectral function upon (in eV)\n\
  [start end steps]. e.g. -15 10 0.5 for an energy step every 0.5 eV between -15 eV and 10 eV. \n')
freq = N.arange(N.float(user_input.split()[0]),N.float(user_input.split()[1]),N.float(user_input.split()[2]))
freq = freq/Ha2eV # From eV to Ha

# Temperature
user_input = raw_input('Enter the temperature for the spectral function A_nk(omega,T) [in K]\n')
temperature = N.float(user_input)

# Band index
user_input = raw_input('Enter the band index for the spectral function A_nk(omega,T)\n')
try:
  band = N.int(user_input)
except ValueError:
  raise Exception('The value you enter is not an integer!')

# Get the nb of random Q-points from user 
user_input = raw_input('Enter the number of random Q-points you have\n')
try:
  nbQ = int(user_input)
except ValueError:
  raise Exception('The value you enter is not an integer!')

# Get the path of the DDB files from user
user_input = raw_input('Enter the name of the %s DDB files separated by a space\n' %nbQ)
if len(user_input.split()) != nbQ:
  raise Exception("You sould provide %s DDB files" %nbQ)
else:
  DDB_files = user_input.split()

# Test if the first file is at the Gamma point
DDBtmp = system(directory='.',filename=DDB_files[0])
if N.allclose(DDBtmp.iqpt,[0.0,0.0,0.0]) == False:
  raise Exception('The first Q-point is not Gamma!')

# Choose a k-point in the list below:
print 'Choose a k-point number in the list below for A_nk(omega,T)\n'
for ii in N.arange(DDBtmp.nkpt):
  print '%s) %s' % (ii,DDBtmp.kpt[ii,:])
user_input = raw_input('Enter the number of the k-point you want to analyse\n')
try:
  kpt = N.int(user_input)
except ValueError:
  raise Exception('The value you enter is not an integer!')

# Get the path of the eigq files from user
user_input = raw_input('Enter the name of the %s eigq files separated by a space\n' %nbQ)
if len(user_input.split()) != nbQ:
  raise Exception("You sould provide %s DDB files" %nbQ)
else:
  eigq_files = user_input.split()

# Get the path of the EIGR2D files from user
user_input = raw_input('Enter the name of the %s EIGR2D files separated by a space\n' %nbQ)
if len(user_input.split()) != nbQ:
  raise Exception("You sould provide %s DDB files" %nbQ)
else:
  EIGR2D_files = user_input.split()

# Get the path of the FAN files from user
user_input = raw_input('Enter the name of the %s FAN files separated by a space\n' %nbQ)
if len(user_input.split()) != nbQ:
  raise Exception("You sould provide %s DDB files" %nbQ)
else:
  FAN_files = user_input.split()

# Take the EIG at Gamma
user_input = raw_input('Enter the name of the unperturbed EIG.nc file at Gamma\n')
if len(user_input.split()) != 1:
  raise Exception("You sould only provide 1 file")
else:
  eig0 = system(directory='.',filename=user_input)

N.arange# Find the degenerate eigenstates
DDB = system(directory='.',filename=DDB_files[0])
degen =  N.zeros((DDB.nkpt,DDB.nband),dtype=int)
for ikpt in N.arange(DDB.nkpt):
  count = 0
  for iband in N.arange(DDB.nband):
    if iband != DDB.nband-1:
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

# Read the EIGR2D file at Gamma and save it in ddw_save
EIGR2D = system()
EIGR2D.__init__(directory='.',filename=EIGR2D_files[0])
ddw_save = N.zeros((3,EIGR2D.natom,3,EIGR2D.natom),dtype=complex)
for iatom1 in N.arange(EIGR2D.natom):
  for iatom2 in N.arange(EIGR2D.natom):
    for idir1 in N.arange(3):
      for idir2 in N.arange(3):
        ddw_save[idir1,iatom1,idir2,iatom2] = EIGR2D.EIG2D[kpt,band-1,idir1,iatom1,idir2,iatom2]

# Create the random Q-integration (wtq=1/nqpt):
wtq = N.ones((nbQ))
wtq = wtq*(1.0/nbQ)

# Parallelize the work over cpus
pool = multiprocessing.Pool(processes=nb_cpus)

total = pool.map(FanDDW, zip(wtq,eigq_files,DDB_files,EIGR2D_files,FAN_files))
FanDDW_corr = sum(total)

# Computation of the self-energy and the spectral function at a given temperature.
print 'The dft eigenenergy is e_nk = ',eig0.EIG[0,kpt,band-1]*Ha2eV
self_energy = N.zeros((len(freq)),dtype=complex)
spectral = N.zeros((len(freq)))
index = 0
with open(output,"w") as O:
  O.write('# Spectral function at T = '+str(temperature)+' for band = '+str(band)+' and kpt = '+str(DDBtmp.kpt[kpt,:])+'\n')
  O.write('# Omega [eV]  Spectral function\n')
  for ifreq in freq:
    self_energy[index] = 1.0/(ifreq-eig0.EIG[0,kpt,band-1]-FanDDW_corr[index])
    spectral[index]= (1.0/N.pi)*N.abs((self_energy[index]).imag)
    O.write(str(ifreq*Ha2eV)+' '+str(spectral[index])+'\n')
    index +=1

# Make a matplotlibplot
#plt.plot(freq*Ha2eV,spectral)
#plt.ylabel('Spectral function')
#plt.xlabel('Energy [eV]')
#plt.show()

# Report wall time
end = datetime.now()
print 'End on %s/%s/%s at %s h %s ' %(end.day,end.month,end.year,end.hour,end.minute)

runtime = end - start
print "Runtime: %s seconds (or %s minutes)" %(runtime.seconds,float(runtime.seconds)/60.0)













