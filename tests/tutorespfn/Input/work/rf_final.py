# Author: S. Ponc\'e + Y. Gillet
# Date: 30/04/2013 -- 11/09/2014 -- 07/08/2015 -- 21/12/2020
# Version 1.5
# Classes needed for the temperature_final.py script
# 2015 : Spin + Lifetime coding
# 2020 : Port to python3
# 2024 : Port to python3.6

#from __future__ import division, print_function
import numpy as N
from numpy import zeros
#from numpy import complex, float
import itertools as Iter
from functools import partial
import multiprocessing
import netCDF4 as nc
import sys
import os

# Variables
tol6 = 1E-6
tol8 = 1E-8
Ha2eV = 27.21138386
kb_HaK = 3.1668154267112283e-06

###########
# CLASSES #
###########

class system:
  natom  = 0
  ntypat = 0
  nkpt   = 0
  kpt    = None
  Kptns  = None
  EIG    = None
  nband  = 0
  acell  = None
  occ    = None
  amu    = None
  rprim  = N.empty((3, 3))
  iqpt   = None
  IFC    = None
  filename     = None
  filefullpath = None
  def __init__(self,directory=None,filename=None):
    if filename == None:return
    if directory == None:directory='.'
    self.filename = filename
    self.filefullpath = '%s/%s' %(directory,filename)
    if self.filefullpath[-4:] == '_DDB':
      self.DDB_file_open(self.filefullpath)
    if self.filefullpath[-10:] == '_EIGR2D.nc' or self.filefullpath[-10:] == '_EIGI2D.nc':
      self.EIG2Dnc_file_open(self.filefullpath)
    if self.filefullpath[-7:] == '_EIG.nc':
      self.EIG_file_open(self.filefullpath)
    if self.filefullpath[-4:] == '_EIG':
      raise Exception('Please provide a netCDF _EIG.nc file!\n\
         This is mandatory for good accuracy.' )
    if self.filefullpath[-7:] == '_GKK.nc':
      self.GKKnc_file_open(self.filefullpath)
    if self.filefullpath[-6:] == '_EP.nc':
      self.EP_file_open(self.filefullpath)

# Read _EP.nc file
  def EP_file_open(self,filefullpath):
    if not (os.path.isfile(filefullpath)):
      raise Exception('The file "%s" does not exists!' %filefullpath)
    root = nc.Dataset(filefullpath,'r')
    self.natom = len(root.dimensions['number_of_atoms'])
    self.nkpt = len(root.dimensions['number_of_kpoints'])
    self.nband = len(root.dimensions['max_number_of_states'])
    self.ntemp = len(root.dimensions['number_of_temperature'])
    self.nsppol = len(root.dimensions['number_of_spins'])
    self.nbQ = len(root.dimensions['number_of_qpoints'])
    self.temp = root.variables['temperature'][:]
    self.occ = root.variables['occupations'][:,:,:] # number_of_spins, number_of_kpoints, max_number_of_states
    self.kpt = root.variables['reduced_coordinates_of_kpoints'][:,:]
    self.eigenvalues = root.variables['eigenvalues'][:,:,:] #number_of_spins, number_of_kpoints, max_number_of_states
    self.rprimd = root.variables['primitive_vectors'][:,:]
    self.zpm = root.variables['zero_point_motion'][:,:,:,:,:] # nsppol, number_of_temperature,
                                                   # number_of_kpoints, max_number_of_states, cplex
    root.close()

# Read _EIG.nc file
  def EIG_file_open(self,filefullpath):
    if not (os.path.isfile(filefullpath)):
      raise Exception('The file "%s" does not exists!' %filefullpath)
    root = nc.Dataset(filefullpath,'r')
    self.EIG = root.variables['Eigenvalues'][:,:,:] # nsppol,nkpt,nband
    self.Kptns = root.variables['Kptns'][:,:]
    NBandK = root.variables['NBandK'][:]
    self.nband =  int(NBandK[0,0])
    root.close()

# Open the Fan.nc file and read it
  def GKKnc_file_open(self,filefullpath):
    if not (os.path.isfile(filefullpath)):
      raise Exception('The file "%s" does not exists!' %filefullpath)
    root = nc.Dataset(filefullpath,'r')
    self.natom = len(root.dimensions['number_of_atoms'])
    self.nkpt = len(root.dimensions['number_of_kpoints'])
    self.nband = len(root.dimensions['max_number_of_states'])
    self.nsppol = len(root.dimensions['number_of_spins'])
    self.occ = root.variables['occupations'][:,:,:] # number_of_spins, number_of_kpoints, max_number_of_states
    GKKtmp = root.variables['second_derivative_eigenenergies_actif'][:,:,:,:,:] #max_number_of_states,number_of_atoms,
                                       # number_of_cartesian_directions, number_of_kpoints, product_mband_nsppol*2
    GKKtmp2 = N.einsum('ijkno->nokji', GKKtmp)
    #GKKtmp2(nkpt,nband*nsppol*2,3,natom,nband)
    GKKtmp3 = GKKtmp2[:, ::2, ...]  # Slice the even numbers
    GKKtmp4 = GKKtmp2[:, 1::2, ...] # Slice the odd numbers
    self.GKK = 1j*GKKtmp4
    self.GKK += GKKtmp3
    self.GKK_bis = N.reshape(self.GKK,(self.nkpt,self.nsppol,self.nband,3,self.natom,self.nband))
    self.eigenvalues = root.variables['eigenvalues'][:,:,:] #number_of_spins, number_of_kpoints, max_number_of_states
    self.kpt = root.variables['reduced_coordinates_of_kpoints'][:,:]
    self.iqpt = root.variables['current_q_point'][:]
    self.wtq = root.variables['current_q_point_weight'][:]
    self.rprimd = root.variables['primitive_vectors'][:,:]
    root.close()

# Open the EIG2D.nc file and read it
  def EIG2Dnc_file_open(self,filefullpath):
    if not (os.path.isfile(filefullpath)):
      raise Exception('The file "%s" does not exists!' %filefullpath)
    root = nc.Dataset(filefullpath,'r')
    self.natom = len(root.dimensions['number_of_atoms'])
    self.nkpt = len(root.dimensions['number_of_kpoints'])
    self.nband = len(root.dimensions['maximum_number_of_bands'])
    self.nsppol = len(root.dimensions['number_of_spins'])
    self.occ = root.variables['occupations'][:,:,:] # number_of_spins, number_of_kpoints, max_number_of_states
    group = root.groups['d2eig']
    EIG2Dtmp = group.variables['matrix_values'][:,:,:,:,:,:,:,:,:]#number_of_d2eig_blocks, number_of_spins, number_of_kpoints, maximum_number_of_bands, number_of_perturbations, number_of_cartesian_directions, number_of_perturbations, number_of_cartesian_directions, cplex
    #EIG2Dtmp = root.variables['second_derivative_eigenenergies'][:,:,:,:,:,:,:] #number_of_atoms,
                                       # number_of_cartesian_directions, number_of_atoms, number_of_cartesian_directions,
                                       # number_of_kpoints, product_mband_nsppol, cplex
    # SP: we assume only 1 q-points (no merge) per file and only 1 spin channel
    EIG2Dtmp2 = EIG2Dtmp[0,0,:,:,:,:,:,:,:] # (nkpt, mband, npert, 3, npert, 3, 2)
    # From (nkpt, mband, npert, 3, npert, 3, 2) --> (nkpt,mband*nsppol,3,natom,3,natom,2)
    EIG2Dtmp3 = N.einsum('ijklmno->ijlknmo', EIG2Dtmp2)

    self.EIG2D = 1j*EIG2Dtmp3[...,1]
    self.EIG2D += EIG2Dtmp3[...,0]
    #EIG2D(nkpt,mband*nsppol,3,natom,3,natom)
    self.EIG2D_bis = N.reshape(self.EIG2D,(self.nkpt,self.nsppol,self.nband,3,self.natom,3,self.natom))
    #EIG2D_bis(nkpt,nband,nsppol,3,natom,3,natom)
    #self.eigenvalues = root.variables['eigenvalues'][:,:,:] #number_of_spins, number_of_kpoints, max_number_of_states
    self.kpt = root.variables['reduced_coordinates_of_kpoints'][:,:]
    self.iqpt = root.groups['d2eig'].variables['reduced_coordinates_of_qpoints'][:]
    #self.iqpt = root.variables['current_q_point'][:]
    #self.wtq = root.variables['current_q_point_weight'][:]
    #self.rprimd = root.variables['primitive_vectors'][:,:]
    root.close()

# Open the DDB file and read it
  def DDB_file_open(self,filefullpath):
    if not (os.path.isfile(filefullpath)):
      raise Exception('The file "%s" does not exists!' %filefullpath)
    with open(filefullpath,'r') as DDB:
      Flag = 0
      Flag2 = False
      Flag3 = False
      ikpt = 0
      typatdone = 0
      for line in DDB:
        if line.find('natom') > -1:
          self.natom = int(line.split()[1])
        if line.find('nkpt') > -1:
          self.nkpt = int(line.split()[1])
          self.kpt  = zeros((self.nkpt,3))
        if line.find('ntypat') > -1:
          self.ntypat = int(line.split()[1])
        if line.find('nband') > -1:
          self.nband = int(line.split()[1])
        if line.find('acell') > -1:
          line = line.replace('D','E')
          tmp = line.split()
          self.acell = [float(tmp[1]),float(tmp[2]),float(tmp[3])]
        if Flag2:
          line = line.replace('D','E')
          for ii in N.arange(3,self.ntypat):
            self.amu[ii] = float(line.split()[ii-3])
            Flag2 = False
        if line.find('amu') > -1:
          line = line.replace('D','E')
          self.amu = zeros((self.ntypat))
          if self.ntypat > 3:
            for ii in N.arange(3):
              self.amu[ii] = float(line.split()[ii+1])
              Flag2 = True
          else:
            for ii in N.arange(self.ntypat):
              self.amu[ii] = float(line.split()[ii+1])
        if line.find(' kpt ') > -1:
          line = line.replace('D','E')
          tmp = line.split()
          self.kpt[0,0:3] = [float(tmp[1]),float(tmp[2]),float(tmp[3])]
          ikpt = 1
          continue
        if ikpt < self.nkpt and ikpt > 0:
          line = line.replace('D','E')
          tmp = line.split()
          self.kpt[ikpt,0:3] = [float(tmp[0]),float(tmp[1]),float(tmp[2])]
          ikpt += 1
          continue
        if Flag == 2:
          line = line.replace('D','E')
          tmp = line.split()
          self.rprim[2,0:3] = [float(tmp[0]),float(tmp[1]),float(tmp[2])]
          Flag = 0
        if Flag == 1:
          line = line.replace('D','E')
          tmp = line.split()
          self.rprim[1,0:3] = [float(tmp[0]),float(tmp[1]),float(tmp[2])]
          Flag = 2
        if line.find('rprim') > -1:
          line = line.replace('D','E')
          tmp = line.split()
          self.rprim[0,0:3] = [float(tmp[1]),float(tmp[2]),float(tmp[3])]
          Flag = 1
        if Flag3:
          line = line.replace('D','E')
          if (self.natom-typatdone)*1.0/12 < 1.001:
            for ii in N.arange(self.natom-typatdone):
              self.typat[typatdone+ii] = float(line.split()[ii])
              Flag3 = False
          else:
            for ii in N.arange(12):
              self.typat[typatdone+ii] = float(line.split()[ii])
            typatdone += 12
        if line.find(' typat') > -1:
          self.typat = zeros((self.natom))
          if self.natom > 12:
            for ii in N.arange(12):
              self.typat[ii] = float(line.split()[ii+1])
              Flag3 = True
              typatdone = 12
          else:
            for ii in N.arange(self.natom):
              self.typat[ii] = float(line.split()[ii+1])
        # Read the actual d2E/dRdR matrix
        if Flag == 3:
          line = line.replace('D','E')
          tmp = line.split()
          self.IFC[int(tmp[0])-1,int(tmp[1])-1,int(tmp[2])-1,int(tmp[3])-1] = \
            complex(float(tmp[4]),float(tmp[5]))
        # Read the current Q-point
        if line.find('qpt') > -1:
          line = line.replace('D','E')
          tmp = line.split()
          self.iqpt = [float(tmp[1]),float(tmp[2]),float(tmp[3])]
          Flag = 3
          self.IFC = zeros((3,self.natom,3,self.natom),dtype=complex)

#################################################
# Usefull definition to avoid code duplications #
#################################################
def compute_dynmat(DDB):
# Retrive the amu for each atom
  amu = zeros(DDB.natom)
  for ii in N.arange(DDB.natom):
    jj = DDB.typat[ii].astype(int)
    amu[ii] = DDB.amu[jj-1]
# Calcul of gprimd from rprimd
  rprimd = DDB.rprim*DDB.acell
  gprimd = N.linalg.inv(N.matrix(rprimd))
# Transform from 2nd-order matrix (non-cartesian coordinates,
# masses not included, asr not included ) from DDB to
# dynamical matrix, in cartesian coordinates, asr not imposed.
  IFC_cart = zeros((3,DDB.natom,3,DDB.natom),dtype=complex)
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
  Dyn_mat = zeros((3*DDB.natom,3*DDB.natom),dtype=complex)
  while ipert1 < 3*DDB.natom:
    for ii in N.arange(DDB.natom):
      for dir1 in N.arange(3):
        ipert2 = 0
        while ipert2 < 3*DDB.natom:
          for jj in N.arange(DDB.natom):
            for dir2 in N.arange(3):
              Dyn_mat[ipert1,ipert2] = IFC_cart[dir1,ii,dir2,jj]*(5.4857990965007152E-4)/ \
                   N.sqrt(amu[ii]*amu[jj])
              ipert2 += 1
        ipert1 += 1
# Hermitianize the dynamical matrix
  dynmat = N.matrix(Dyn_mat)
  dynmat = 0.5*(dynmat + dynmat.transpose().conjugate())

# Solve the eigenvalue problem with linear algebra (Diagonalize the matrix)
  [eigval,eigvect]=N.linalg.eigh(Dyn_mat)

# Orthonormality relation
  ipert = 0
  for ii in N.arange(DDB.natom):
    for dir1 in N.arange(3):
     eigvect[ipert] = (eigvect[ipert])*N.sqrt(5.4857990965007152E-4/amu[ii])
     ipert += 1
  kk = 0
  for jj in eigval:
    if jj < 0.0:
      print("WARNING: An eigenvalue is negative with value: ",jj," ... but proceed with value 0.0")
      eigval[kk] = 0.0
      kk += 1
    else:
      kk += 1
  omega = N.sqrt(eigval) #*5.4857990965007152E-4)
#  print "omega",omega

# The acoustic phonon at Gamma should NOT contribute because they should be zero.
# Moreover with the translational invariance the ZPM will be 0 anyway for these
# modes but the FAN and DDW will have a non physical value. We should therefore
# neglect these values.
#  if N.allclose(DDB.iqpt,[0.0,0.0,0.0]) == True:
#    omega[0] = 0.0
#    omega[1] = 0.0
#    omega[2] = 0.0

  return omega,eigvect,gprimd

# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------

def get_reduced_displ(natom,eigvect,omega,gprimd):
  displ_FAN =  zeros((3,3),dtype=complex)
  displ_DDW =  zeros((3,3),dtype=complex)
  displ_red_FAN2 = zeros((3*natom,natom,natom,3,3),dtype=complex)
  displ_red_DDW2 = zeros((3*natom,natom,natom,3,3),dtype=complex)
  for imode in N.arange(3*natom): #Loop on perturbation (6 for 2 atoms)
    if omega[imode].real > tol6:
      for iatom1 in N.arange(natom):
        for iatom2 in N.arange(natom):
          for idir1 in N.arange(0,3):
            for idir2 in N.arange(0,3):
              displ_FAN[idir1,idir2] = eigvect[3*iatom2+idir2,imode].conj()\
                 *eigvect[3*iatom1+idir1,imode]/(2.0*omega[imode].real)
              displ_DDW[idir1,idir2] = (eigvect[3*iatom2+idir2,imode].conj()\
                 *eigvect[3*iatom2+idir1,imode]+eigvect[3*iatom1+idir2,imode].conj()\
                 *eigvect[3*iatom1+idir1,imode])/(4.0*omega[imode].real)
              # Now switch to reduced coordinates in 2 steps (more efficient)
          tmp_displ_FAN = zeros((3,3),dtype=complex)
          tmp_displ_DDW = zeros((3,3),dtype=complex)
          for idir1 in N.arange(3):
            for idir2 in N.arange(3):
              tmp_displ_FAN[:,idir1] = tmp_displ_FAN[:,idir1]+displ_FAN[:,idir2]*gprimd[idir2,idir1]
              tmp_displ_DDW[:,idir1] = tmp_displ_DDW[:,idir1]+displ_DDW[:,idir2]*gprimd[idir2,idir1]
          displ_red_FAN = zeros((3,3),dtype=complex)
          displ_red_DDW = zeros((3,3),dtype=complex)
          for idir1 in N.arange(3):
            for idir2 in N.arange(3):
              displ_red_FAN[idir1,:] = displ_red_FAN[idir1,:] + tmp_displ_FAN[idir2,:]*gprimd[idir2,idir1]
              displ_red_DDW[idir1,:] = displ_red_DDW[idir1,:] + tmp_displ_DDW[idir2,:]*gprimd[idir2,idir1]

          displ_red_FAN2[imode,iatom1,iatom2,:,:] = displ_red_FAN[:,:]
          displ_red_DDW2[imode,iatom1,iatom2,:,:] = displ_red_DDW[:,:]
  return displ_red_FAN2,displ_red_DDW2

# ----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------

def make_average(nsppol,nkpt,nband,degen,total_corr,temp=False):
  if temp:
    for ispin in N.arange(nsppol):
      for ikpt in N.arange(nkpt):
        count = 0
        iband = 0
        while iband < nband:
          if iband < nband-2:
            if ((degen[ispin,ikpt,iband] == degen[ispin,ikpt,iband+1]) and (degen[ispin,ikpt,iband] == degen[ispin,ikpt,iband+2])):
              total_corr[:,:,ispin,ikpt,iband] = (total_corr[:,:,ispin,ikpt,iband]+total_corr[:,:,ispin,ikpt,iband+1]+total_corr[:,:,ispin,ikpt,iband+2])/3
              total_corr[:,:,ispin,ikpt,iband+1] = total_corr[:,:,ispin,ikpt,iband]
              total_corr[:,:,ispin,ikpt,iband+2] = total_corr[:,:,ispin,ikpt,iband]
              iband += 3
              continue
          if iband <  nband-1:
            if (degen[ispin,ikpt,iband] == degen[ispin,ikpt,iband+1]):
              total_corr[:,:,ispin,ikpt,iband] = (total_corr[:,:,ispin,ikpt,iband]+total_corr[:,:,ispin,ikpt,iband+1])/2
              total_corr[:,:,ispin,ikpt,iband+1]=total_corr[:,:,ispin,ikpt,iband]
              iband +=2
              continue
          iband += 1
  else:
    for ispin in N.arange(nsppol):
      for ikpt in N.arange(nkpt):
        count = 0
        iband = 0
        while iband < nband:
          if iband < nband-2:
            if ((degen[ispin,ikpt,iband] == degen[ispin,ikpt,iband+1]) and (degen[ispin,ikpt,iband] == degen[ispin,ikpt,iband+2])):
              total_corr[:,ispin,ikpt,iband] = (total_corr[:,ispin,ikpt,iband]+total_corr[:,ispin,ikpt,iband+1]+\
                                       total_corr[:,ispin,ikpt,iband+2])/3
              total_corr[:,ispin,ikpt,iband+1] = total_corr[:,ispin,ikpt,iband]
              total_corr[:,ispin,ikpt,iband+2] = total_corr[:,ispin,ikpt,iband]
              iband += 3
              continue
          if iband <  nband-1:
            if (degen[ispin,ikpt,iband] == degen[ispin,ikpt,iband+1]):
              total_corr[:,ispin,ikpt,iband] = (total_corr[:,ispin,ikpt,iband]+total_corr[:,ispin,ikpt,iband+1])/2
              total_corr[:,ispin,ikpt,iband+1]=total_corr[:,ispin,ikpt,iband]
              iband +=2
              continue
          iband += 1
  return total_corr


# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------

def get_bose(natom,omega,temp_info):
  bose = N.array(zeros((3*natom,len(temp_info))))
  for imode in N.arange(3*natom): #Loop on perturbation (6 for 2 atoms)
    if omega[imode].real > tol6:
      tt = 0
      for T in temp_info:
        if T < tol6:
          bose[imode,tt] = 0.0
        else:
          bose[imode,tt] = 1.0/(N.exp(omega[imode].real/(kb_HaK*T))-1)
        tt += 1
  #print bose[:,0]
  return bose

# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------

def gaussian(ifreq,omegatmp,gaussian_smearing):
  x = (ifreq - omegatmp)/gaussian_smearing
  gaussian = 1.0/(N.sqrt(N.pi)) * (1.0/gaussian_smearing) * N.exp(-(x**2))
  return gaussian

def fermidirac(freq,center,smearing):
  x = (freq-center)/smearing
  value = 0.25/smearing/(N.cosh(x/2.0))**2
  return value

def lorentzian(freq,center,smearing):
  x = (freq-center)/smearing
  value = 1.0/(N.pi)/smearing/(x**2+1)
  return value


#####################
# Compute temp. dep #
#####################

# Compute the dynamical ZPR with temperature dependence
def dynamic_zpm_temp(arguments,ddw_save,ddw_save2,type,temp_info,smearing,eig0,degen,energy,gaussian_smearing):

  if type == 1 or type == 2:
    nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files,GKK_files = arguments
    GKKterm = system(directory='.',filename=GKK_files)
    GKK = GKKterm.GKK_bis
  elif type == 3:
    nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files = arguments

  DDB = system(directory='.',filename=DDB_files)
  EIGR2D = system(directory='.',filename=EIGR2D_files)
  ntemp = len(temp_info)
  if type == 1 or type == 2:
    total_corr = zeros((4+2*len(energy),ntemp,EIGR2D.nsppol,EIGR2D.nkpt,EIGR2D.nband),dtype=complex)
  elif type == 3:
    total_corr = zeros((4+len(energy),ntemp,EIGR2D.nsppol,EIGR2D.nkpt,EIGR2D.nband),dtype=complex)
  eigq = system(directory='.',filename=eigq_files)

# If the calculation is on a Homogenous q-point mesh
# retreve the weight of the q-point
  if (wtq == 0):
    wtq = EIGR2D.wtq
    wtq = wtq[0]

# Current Q-point calculated
  print("Q-point: ",nbqpt," with wtq =",wtq," and reduced coord.",EIGR2D.iqpt)
  current = multiprocessing.current_process()
  file_name = str('PYLOG_')+str(current.pid)
  if os.path.isfile(file_name) :
    with open(file_name,'a') as F:
      F.write("Q-point: "+str(nbqpt)+" with wtq ="+str(wtq)+" and reduced coord."+str(EIGR2D.iqpt)+"\n")
  else:
    with open(file_name,'w') as F:
      F.write("Q-point: "+str(nbqpt)+" with wtq ="+str(wtq)+" and reduced coord."+str(EIGR2D.iqpt)+"\n")


# Find phonon freq and eigendisplacement from _DDB
  omega,eigvect,gprimd=compute_dynmat(DDB)

# Compute the displacement = eigenvectors of the DDB.
# Due to metric problem in reduce coordinate we have to work in cartesian
# but then go back to reduce because our EIGR2D matrix elements are in reduced coord.
  fan_corr =  zeros((ntemp,EIGR2D.nkpt,EIGR2D.nband,EIGR2D.nsppol),dtype=complex)
  ddw_corr = zeros((ntemp,EIGR2D.nkpt,EIGR2D.nband,EIGR2D.nsppol),dtype=complex)
  fan_add = N.array(zeros((ntemp,EIGR2D.nkpt,EIGR2D.nband,EIGR2D.nsppol),dtype=complex))
  ddw_add = N.array(zeros((ntemp,EIGR2D.nkpt,EIGR2D.nband,EIGR2D.nsppol),dtype=complex))

  bose = get_bose(EIGR2D.natom,omega,temp_info)
  # bose (imode, itemp)

# Get reduced displacement (scaled with frequency)
  displ_red_FAN2,displ_red_DDW2 = get_reduced_displ(EIGR2D.natom,eigvect,omega,gprimd)
  # displ_red(mode,atom1,atom2,dir1,dir2)

# Einstein sum make the vector matrix multiplication ont the correct indices
  fan_corrQ = N.einsum('iojklmn,plnkm->pijo',EIGR2D.EIG2D_bis,displ_red_FAN2)
  ddw_corrQ = N.einsum('iojklmn,plnkm->pijo',ddw_save,displ_red_DDW2)
  # fan_corrQ(mode,kpt,band,spin)

  # Sum over the modes with bose + reshape to (itemp,ispin,ikpt,iband)
  fan_corr = N.einsum('ijkl,im->mljk',fan_corrQ,2*bose+1.0)
  ddw_corr = N.einsum('ijkl,im->mljk',ddw_corrQ,2*bose+1.0)

  omegatmp = omega[:].real # imode
  if type == 3:
    imag_fan_add = 0.0
  elif type == 1 or type == 2:
    print("Now compute active space ...")

#   Now compute active space
    # sum over atom1,dir1
    temp = N.einsum('iqjklm,nlokp->ijmnpoq',GKK,displ_red_FAN2)
    fan_addQ = N.einsum('ijklmnq,iqjmnk->ijklq',temp,N.conjugate(GKK))
    # fan_addQ(nkpt,nband,nband,imode,ispin)
    temp = N.einsum('iqjklm,nlokp->ijmnpoq',ddw_save2,displ_red_DDW2)
    ddw_addQ = N.einsum('ijklmnq,iqjmnk->ijklq',temp,N.conjugate(ddw_save2))
    # fan_addQ(nkpt,nband,nband,imode,ispin)

    occtmp = EIGR2D.nsppol*EIGR2D.occ[:,:,:]/2 # jband # should be 1 !
    delta_E_ddw = N.einsum('lij,k->lijk',eig0[:,:,:].real,N.ones(EIGR2D.nband)) - \
              N.einsum('lij,k->likj',eig0[:,:,:].real,N.ones(EIGR2D.nband)) - \
              N.einsum('i,ljk->ljik',N.ones((EIGR2D.nband)),(2*occtmp-1))*smearing*1j # spin,ikpt,iband,jband

    ddw_tmp = N.einsum('ijkln,lm->mijkn',ddw_addQ,2*bose+1.0) # itemp,ikpt,iband,jband,ispin
    ddw_add = N.einsum('ijklm,mjkl->imjk',ddw_tmp,1.0/delta_E_ddw) # temp,spin,ikpt,iband
    delta_E = N.einsum('lij,k->lijk',eig0[:,:,:].real,N.ones(EIGR2D.nband)) - \
              N.einsum('lij,k->likj',eigq.EIG[:,:,:].real,N.ones(EIGR2D.nband)) # spin,ikpt,iband,jband
    delta_E_sm = N.einsum('i,ljk->ljik',N.ones((EIGR2D.nband)),(2*occtmp-1))*smearing*1j # spin,ikpt,iband,jband
    num1 = N.einsum('ij,mkl->mkijl',bose,N.ones((EIGR2D.nsppol,EIGR2D.nkpt,EIGR2D.nband))) +1.0 \
          - N.einsum('ij,mkl->mkijl',N.ones((3*EIGR2D.natom,ntemp)),occtmp) # spin,k,mod,temp,band # bef was (imode,tmp,band)
    deno1 = N.einsum('mijk,l->mijkl',delta_E,N.ones(3*EIGR2D.natom),dtype=complex)

    if type==1: # dynamic
      deno1 -= N.einsum('mijk,l->mijkl',N.ones((EIGR2D.nsppol,EIGR2D.nkpt,EIGR2D.nband,EIGR2D.nband)),omegatmp) #spin,ikpt,iband,jband,imode

    imag_part1 = N.pi*gaussian(deno1,0.0,gaussian_smearing)
    #imag_part1 = N.pi*fermidirac(deno1,0.0,gaussian_smearing)
    #imag_part1 = N.pi*lorentzian(deno1,0.0,gaussian_smearing)
    #imag_part1 = gaussian_smearing/(deno1*deno1 + gaussian_smearing*gaussian_smearing)
    deno1 += N.einsum('mijk,l->mijkl',delta_E_sm,N.ones(3*EIGR2D.natom),dtype=complex)

    div1 = N.einsum('ijklm,ijnmk->iklmjn',num1,1.0/deno1) # (spin,k,mod,temp,jband)/(spin,ikpt,iband,jband,mode) => (ispin,imod,tmp,jband,ikpt,iband)

    num2 = N.einsum('ij,mkl->mkijl',bose,N.ones((EIGR2D.nsppol,EIGR2D.nkpt,EIGR2D.nband))) \
          + N.einsum('ij,mkl->mkijl',N.ones((3*EIGR2D.natom,ntemp)),occtmp) #imode,tmp,jband
    deno2 = N.einsum('mijk,l->mijkl',delta_E,N.ones((3*EIGR2D.natom),dtype=complex))

    if type==1: # dynamic
      deno2 += N.einsum('mijk,l->mijkl',N.ones((EIGR2D.nsppol,EIGR2D.nkpt,EIGR2D.nband,EIGR2D.nband)),omegatmp) #spin,ikpt,iband,jband,imode

    imag_part2 = N.pi*gaussian(deno2,0.0,gaussian_smearing)
    #imag_part2 = N.pi*fermidirac(deno2,0.0,gaussian_smearing)
    #imag_part2 = N.pi*lorentzian(deno2,0.0,gaussian_smearing)
    #imag_part2 = gaussian_smearing/(deno2*deno2 + gaussian_smearing*gaussian_smearing)
    deno2 -= N.einsum('mijk,l->mijkl',delta_E_sm,N.ones(3*EIGR2D.natom),dtype=complex)

    div2 = N.einsum('ijklm,ijnmk->iklmjn',num2,1.0/deno2) # (spin,k,mod,temp,jband)/(spin,ikpt,iband,jband,mode) => (ispin,imod,tmp,jband,ikpt,iband)

    fan_add = N.einsum('ijklq,qlmkij->mqij',fan_addQ,div1+div2) # (k,iband,jband,imod,ispin) * (spin,imod,tmp,jband,ikpt,iband) => (temp,ispin,ikpt,iband)


    imag_div1 = N.einsum('ijklm,ijnmk->iklmjn',num1,imag_part1) # (spin,k,mod,temp,jband)/(spin,ikpt,iband,jband,mode) => (ispin,imod,tmp,jband,ikpt,iband)
    imag_div2 = N.einsum('ijklm,ijnmk->iklmjn',num2,imag_part2) # (spin,k,mod,temp,jband)/(spin,ikpt,iband,jband,mode) => (ispin,imod,tmp,jband,ikpt,iband)
    imag_fan_add = N.einsum('ijklq,qlmkij->mqij',fan_addQ,imag_div1+imag_div2) # (k,iband,jband,imod,ispin) * (spin,imod,tmp,jband,ikpt,iband) => (temp,ispin,ikpt,iband)

    print("Now compute generalized g2F Eliashberg electron-phonon spectral function ...")

    fan_tmp = N.einsum('ijklm->mijl',fan_addQ) # (ispin,ikpt, iband, imode)
    ddw_tmp = N.einsum('ijklm->mijl',ddw_addQ) # (ispin,ikpt, iband, imode)
    g_kk = fan_tmp - ddw_tmp

    # Eliashberg function
    a2F =  zeros((len(energy),EIGR2D.nsppol,EIGR2D.nkpt,EIGR2D.nband),dtype=complex)
    for ifreq,freq in enumerate(energy):
      for imode in N.arange(3*EIGR2D.natom):
        a2F[ifreq,:,:,:] += g_kk[:,:,:,imode]*gaussian(freq,omegatmp[imode],gaussian_smearing)
    fan_corr += fan_add
    ddw_corr += ddw_add

  # PDOS
  PDOS = zeros((len(energy)),dtype=complex)
  ii = 0
  for ifreq,freq in enumerate(energy):
    for imode in N.arange(3*EIGR2D.natom):
      PDOS[ifreq] += gaussian(freq,omegatmp[imode],gaussian_smearing)

  # From the equations ddw_corr has no physical imaginary part
  ddw_corr.imag = 0.0

  eigen_corr = (fan_corr[:,:,:,:] - ddw_corr[:,:,:,:])*wtq

  total_corr[0,:,:,:,:] = eigen_corr[:,:,:,:]
  total_corr[1,:,:,:,:] = fan_corr[:,:,:,:]*wtq
  total_corr[2,:,:,:,:] = ddw_corr[:,:,:,:]*wtq
  total_corr[3,:,:,:,:] = imag_fan_add*wtq
  total_corr[4:len(energy)+4,0,0,0,0] = PDOS

  if type == 1 or type == 2:
    total_corr[4+len(energy):len(energy)*2+4,0,:,:,:] = a2F

  total_corr = make_average(EIGR2D.nsppol,EIGR2D.nkpt,EIGR2D.nband,degen,total_corr,temp=True)

  return total_corr

#########################################################################################################
# Compute total weigth
def compute_wtq(arguments,type):
  if type == 1 or type == 2:
    nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files,GKK_files = arguments
  elif type == 3:
    nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files = arguments

  EIGR2D = system(directory='.',filename=EIGR2D_files)
# If the calculation is on a Homogenous q-point mesh
# retreve the weight of the q-point
  if (wtq == 0):
    wtq = EIGR2D.wtq
    wtq = wtq[0]
  return wtq


class zpm:
  total_corr = None
  def __init__(self,arguments,ddw_save,ddw_save2,nb_cpus,type,temp_info,smearing,eig0,degen,energy,gaussian_smearing):
    # Parallelize the work over cpus
    pool = multiprocessing.Pool(processes=nb_cpus)
    partial_compute_wtq = partial(compute_wtq,type=type)
    total = pool.imap(partial_compute_wtq,arguments)
    self.total_wtq = sum(total)

    if (type == 1 or type == 2 or type == 3):
      partial_dynamic_zpm_temp =  partial(dynamic_zpm_temp,ddw_save=ddw_save,ddw_save2=ddw_save2,type=type,temp_info=temp_info,\
                                   smearing=smearing,eig0=eig0,degen=degen,energy=energy,gaussian_smearing=gaussian_smearing)
      total = pool.imap(partial_dynamic_zpm_temp,arguments)
      self.total_corr = sum(total)

    pool.close()
    pool.join()

