# Author: Samuel Ponc\'e
# Date: 30/04/2013 -- 11/09/2014
# Version 1.3
# Classes needed for the temperature.py script
# Last devel info: Dynamical coding done

import numpy as N
from numpy import zeros
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
  natom = None
  ntypat = None
  nkpt = None
  kpt = None
  Kptns = None
  EIG = None
  nband = None
  acell = None
  occ = None
  amu = None
  rprim = N.empty((3,3))
  iqpt = None
  IFC = None
  filename = None
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
    self.EIG = root.variables['Eigenvalues'][:,:] 
    self.Kptns = root.variables['Kptns'][:,:]
    NBandK = root.variables['NBandK'][:]
    self.nband =  N.int(NBandK[0,0])
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
    GKKtmp = root.variables['second_derivative_eigenenergies_actif'][:,:,:,:,:] #product_mband_nsppol,number_of_atoms, 
                                       # number_of_cartesian_directions, number_of_kpoints, product_mband_nsppol*2
    GKKtmp2 = zeros((self.nkpt,2*self.nband,3,self.natom,self.nband))
    GKKtmp2 = N.einsum('ijkno->nokji', GKKtmp)
    GKKtmp3 = GKKtmp2[:, ::2, ...]  # Slice the even numbers
    GKKtmp4 = GKKtmp2[:, 1::2, ...] # Slice the odd numbers
    self.GKK = 1j*GKKtmp4
    self.GKK += GKKtmp3
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
    self.nband = len(root.dimensions['max_number_of_states'])
    self.nsppol = len(root.dimensions['number_of_spins'])
    self.occ = root.variables['occupations'][:,:,:] # number_of_spins, number_of_kpoints, max_number_of_states
    EIG2Dtmp = root.variables['second_derivative_eigenenergies'][:,:,:,:,:,:,:] #number_of_atoms, 
                                       # number_of_cartesian_directions, number_of_atoms, number_of_cartesian_directions,
                                       # number_of_kpoints, product_mband_nsppol, cplex
    EIG2Dtmp2 = zeros((self.nkpt,2*self.nband,3,self.natom,3,self.natom,2))
    EIG2Dtmp2 = N.einsum('ijklmno->mnlkjio', EIG2Dtmp)
    self.EIG2D = 1j*EIG2Dtmp2[...,1]
    self.EIG2D += EIG2Dtmp2[...,0]
    self.eigenvalues = root.variables['eigenvalues'][:,:,:] #number_of_spins, number_of_kpoints, max_number_of_states   
    self.kpt = root.variables['reduced_coordinates_of_kpoints'][:,:]
    self.iqpt = root.variables['current_q_point'][:]
    self.wtq = root.variables['current_q_point_weight'][:]
    self.rprimd = root.variables['primitive_vectors'][:,:]
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
      for line in DDB:
        if line.find('natom') > -1:
          self.natom = N.int(line.split()[1])
        if line.find('nkpt') > -1:
          self.nkpt = N.int(line.split()[1])
          self.kpt  = zeros((self.nkpt,3))
        if line.find('ntypat') > -1:
          self.ntypat = N.int(line.split()[1])
        if line.find('nband') > -1:
          self.nband = N.int(line.split()[1])
        if line.find('acell') > -1:
          line = line.replace('D','E')
          tmp = line.split()
          self.acell = [N.float(tmp[1]),N.float(tmp[2]),N.float(tmp[3])]
        if Flag2:
          line = line.replace('D','E')
          for ii in N.arange(3,self.ntypat):
            self.amu[ii] = N.float(line.split()[ii-3])
            Flag2 = False
        if line.find('amu') > -1:
          line = line.replace('D','E')
          self.amu = zeros((self.ntypat))
          if self.ntypat > 3:
            for ii in N.arange(3):
              self.amu[ii] = N.float(line.split()[ii+1])
              Flag2 = True 
          else:
            for ii in N.arange(self.ntypat):
              self.amu[ii] = N.float(line.split()[ii+1])
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
          for ii in N.arange(12,self.natom): 
            self.typat[ii] = N.float(line.split()[ii-12]) 
          Flag3 = False 
        if line.find(' typat') > -1:
          self.typat = zeros((self.natom))
          if self.natom > 12:
            for ii in N.arange(12):
              self.typat[ii] = N.float(line.split()[ii+1])
              Flag3 = True
          else:
            for ii in N.arange(self.natom):
              self.typat[ii] = N.float(line.split()[ii+1])
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
          self.iqpt = [N.float(tmp[1]),N.float(tmp[2]),N.float(tmp[3])]
          Flag = 3
          self.IFC = zeros((3,self.natom,3,self.natom),dtype=complex)

#################################################
# Usefull definition to avoid code duplications #
#################################################
def compute_dynmat(DDB):
# Retrive the amu for each atom
  amu = zeros(DDB.natom)
  for ii in N.arange(DDB.natom):
    jj = DDB.typat[ii]
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
      print "WARNING: An eigenvalue is negative with value: ",jj," ... but proceed with value 0.0"
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

def make_average(nkpt,nband,degen,total_corr,temp=False):
  if temp:
    for ikpt in N.arange(nkpt):
      count = 0
      iband = 0
      while iband < nband:
        if iband < nband-2:
          if ((degen[ikpt,iband] == degen[ikpt,iband+1]) and (degen[ikpt,iband] == degen[ikpt,iband+2])):
            total_corr[:,:,ikpt,iband] = (total_corr[:,:,ikpt,iband]+total_corr[:,:,ikpt,iband+1]+total_corr[:,:,ikpt,iband+2])/3
            total_corr[:,:,ikpt,iband+1] = total_corr[:,:,ikpt,iband]
            total_corr[:,:,ikpt,iband+2] = total_corr[:,:,ikpt,iband]
            iband += 3
            continue
        if iband <  nband-1:
          if (degen[ikpt,iband] == degen[ikpt,iband+1]):
            total_corr[:,:,ikpt,iband] = (total_corr[:,:,ikpt,iband]+total_corr[:,:,ikpt,iband+1])/2
            total_corr[:,:,ikpt,iband+1]=total_corr[:,:,ikpt,iband]
            iband +=2
            continue
        iband += 1
  else:
    for ikpt in N.arange(nkpt):
      count = 0
      iband = 0
      while iband < nband:
        if iband < nband-2:
          if ((degen[ikpt,iband] == degen[ikpt,iband+1]) and (degen[ikpt,iband] == degen[ikpt,iband+2])):
            total_corr[:,ikpt,iband] = (total_corr[:,ikpt,iband]+total_corr[:,ikpt,iband+1]+\
                                     total_corr[:,ikpt,iband+2])/3
            total_corr[:,ikpt,iband+1] = total_corr[:,ikpt,iband]
            total_corr[:,ikpt,iband+2] = total_corr[:,ikpt,iband]
            iband += 3
            continue
        if iband <  nband-1:
          if (degen[ikpt,iband] == degen[ikpt,iband+1]):
            total_corr[:,ikpt,iband] = (total_corr[:,ikpt,iband]+total_corr[:,ikpt,iband+1])/2
            total_corr[:,ikpt,iband+1]=total_corr[:,ikpt,iband]
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




#####################
# Compute temp. dep #
#####################

# Compute the dynamical ZPR only
def dynamic_zpm(arguments,ddw_save,ddw_save2,type,smearing,eig0,degen,energy,gaussian_smearing):
  nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files,GKK_files = arguments
  GKKterm = system(directory='.',filename=GKK_files)
  GKK = GKKterm.GKK
  DDB = system(directory='.',filename=DDB_files)
  EIGR2D = system(directory='.',filename=EIGR2D_files)
  total_corr = zeros((5+2*len(energy),EIGR2D.nkpt,EIGR2D.nband),dtype=complex)
  eigq = system(directory='.',filename=eigq_files)

# If the calculation is on a Homogenous q-point mesh
# retreve the weight of the q-point
  if (wtq == 0):
    wtq = EIGR2D.wtq
    wtq = wtq[0]

# Current Q-point calculated
  print "Q-point: ",nbqpt," with wtq =",wtq," and reduced coord.",EIGR2D.iqpt
  current = multiprocessing.current_process()
#  file_name = str('PYLOG_')+str(current.pid)
#  if os.path.isfile(file_name) :
#    with open(file_name,'a') as F:
#      F.write("Q-point: "+str(nbqpt)+" with wtq ="+str(wtq)+" and reduced coord."+str(EIGR2D.iqpt)+"\n")
#  else:
#    with open(file_name,'w') as F:
#      F.write("Q-point: "+str(nbqpt)+" with wtq ="+str(wtq)+" and reduced coord."+str(EIGR2D.iqpt)+"\n")


# Find phonon freq and eigendisplacement from _DDB
  omega,eigvect,gprimd=compute_dynmat(DDB)

  fan_corr =  zeros((EIGR2D.nkpt,EIGR2D.nband),dtype=complex)
  ddw_corr = zeros((EIGR2D.nkpt,EIGR2D.nband),dtype=complex)
  fan_add = zeros((EIGR2D.nkpt,EIGR2D.nband),dtype=complex)
  ddw_add = zeros((EIGR2D.nkpt,EIGR2D.nband),dtype=complex)

# Get reduced displacement (scaled with frequency)
  displ_red_FAN2,displ_red_DDW2 = get_reduced_displ(EIGR2D.natom,eigvect,omega,gprimd)
# Einstein sum make the vector matrix multiplication ont the correct indices
# fan_corrQ and ddw_corrQ contains the ZPR on stern space.
  fan_corrQ = N.einsum('ijklmn,olnkm->oij',EIGR2D.EIG2D,displ_red_FAN2)
  ddw_corrQ = N.einsum('ijklmn,olnkm->oij',ddw_save,displ_red_DDW2)

  fan_corr = N.sum(fan_corrQ,axis=0)
  ddw_corr = N.sum(ddw_corrQ,axis=0)

  print "Now compute active space ..."

# Now computa active space

#  FAN = N.einsum('ijklm,ijnom->ijklnom',GKK,N.conjugate(GKK))
#  fan_addQ = N.einsum('ijklmno,plnkm->ijop',FAN,displ_red_FAN2) 
  

  temp = N.einsum('ijklm,nlokp->ijmnpo',GKK,displ_red_FAN2)
  fan_addQ = N.einsum('ijklmn,ijmnk->ijkl',temp,N.conjugate(GKK)) 
  temp = N.einsum('ijklm,nlokp->ijmnpo',ddw_save2,displ_red_DDW2)
  ddw_addQ = N.einsum('ijklmn,ijmnk->ijkl',temp,N.conjugate(ddw_save2))

  if type == 1:
    ddw_tmp = N.sum(ddw_addQ,axis=3)
    occtmp = EIGR2D.occ[0,0,:]/2 # jband  
    delta_E = N.einsum('ij,k->ijk',eig0[0,:,:].real,N.ones(EIGR2D.nband)) - \
              N.einsum('ij,k->ikj',eigq.EIG[0,:,:].real,N.ones(EIGR2D.nband)) - \
              N.einsum('ij,k->ijk',N.ones((EIGR2D.nkpt,EIGR2D.nband)),(2*occtmp-1))*smearing*1j # ikpt,iband,jband
    delta_E_ddw = N.einsum('ij,k->ijk',eig0[0,:,:].real,N.ones(EIGR2D.nband)) - \
              N.einsum('ij,k->ikj',eig0[0,:,:].real,N.ones(EIGR2D.nband)) - \
              N.einsum('ij,k->ijk',N.ones((EIGR2D.nkpt,EIGR2D.nband)),(2*occtmp-1))*smearing*1j
    ddw_add = N.einsum('ijk,ijk->ij',ddw_tmp,1.0/delta_E_ddw)
    omegatmp = omega[:].real # imode
    num1 = 1.0-occtmp # jband
    deno1 = N.einsum('ijk,l->ijkl',delta_E,N.ones(3*EIGR2D.natom)) \
          - N.einsum('ijk,l->ijkl',N.ones((EIGR2D.nkpt,EIGR2D.nband,EIGR2D.nband)),omegatmp) #ikpt,iband,jband,imode
    div1 = N.einsum('i,jkil->lijk',num1,1.0/deno1) # (jband)/(ikpt,iband,jband,imode) ==> imode,jband,ikpt,iband
    deno2 = N.einsum('ijk,l->ijkl',delta_E,N.ones(3*EIGR2D.natom)) \
          + N.einsum('ijk,l->ijkl',N.ones((EIGR2D.nkpt,EIGR2D.nband,EIGR2D.nband)),omegatmp) #ikpt,iband,jband,imode
    div2 = N.einsum('i,jkil->lijk',occtmp,1.0/deno2) # (jband)/(ikpt,iband,jband,imode) ==> imode,jband,ikpt,iband
    fan_add = N.einsum('ijkl,lkij->ij',fan_addQ,div1+div2) # ikpt,iband,jband,imod

    print "Now compute generalized g2F Eliashberg electron-phonon spectral function ..." 
 
    fan_tmp = N.sum(fan_addQ,axis=2) # (ikpt, iband, imode)
    ddw_tmp = N.sum(ddw_addQ,axis=2) # (ikpt, iband, imode)
    g_kk = fan_tmp - ddw_tmp
    
    # Generalized Eliashberg function with Gaussian broadening
    a2F =  zeros((len(energy),EIGR2D.nkpt,EIGR2D.nband),dtype=complex)
    PDOS = zeros((len(energy)),dtype=complex)
    ii = 0 
    for ifreq in energy:
      for imode in N.arange(3*EIGR2D.natom):
        a2F[ii,:,:] += g_kk[:,:,imode]*gaussian(ifreq,omegatmp[imode],gaussian_smearing)
        PDOS[ii] += gaussian(ifreq,omegatmp[imode],gaussian_smearing)
      ii += 1
#    a2F_av = make_average(EIGR2D.nkpt,EIGR2D.nband,degen,a2F)

# Correction from active space 
  fan_corr += fan_add
  ddw_corr += ddw_add
  eigen_corr = (fan_corr[:,:] - ddw_corr[:,:])*wtq
# Index 0 contains total ZPR
  total_corr[0,:,:] = eigen_corr[:,:]
# Index 1 contains the FAN correction
  total_corr[1,:,:] = fan_corr[:,:]*wtq
# Index 2 contains the DDW correction
  total_corr[2,:,:] = ddw_corr[:,:]*wtq
# Index 3-energy contains the Eliashberg functions
  total_corr[3:len(energy)+3,:,:] = a2F*wtq 

  total_corr = make_average(EIGR2D.nkpt,EIGR2D.nband,degen,total_corr)

# Index energy+4 contains the PDOS
  total_corr[len(energy)+4:len(energy)*2+4,0,0] = PDOS*wtq

  return total_corr

#########################################################################################################

# Compute the dynamical ZPR with temperature dependence
def dynamic_zpm_temp(arguments,ddw_save,ddw_save2,type,temp_info,smearing,eig0,degen,energy,gaussian_smearing):
  nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files,GKK_files = arguments
  GKKterm = system(directory='.',filename=GKK_files)
  GKK = GKKterm.GKK
  DDB = system(directory='.',filename=DDB_files)
  EIGR2D = system(directory='.',filename=EIGR2D_files)
  total_corr = zeros((5+2*len(energy),len(temp_info),EIGR2D.nkpt,EIGR2D.nband),dtype=complex)
  eigq = system(directory='.',filename=eigq_files)

# If the calculation is on a Homogenous q-point mesh
# retreve the weight of the q-point
  if (wtq == 0):
    wtq = EIGR2D.wtq
    wtq = wtq[0]

# Current Q-point calculated
  print "Q-point: ",nbqpt," with wtq =",wtq," and reduced coord.",EIGR2D.iqpt
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
  fan_corr =  zeros((len(temp_info),EIGR2D.nkpt,EIGR2D.nband),dtype=complex)
  ddw_corr = zeros((len(temp_info),EIGR2D.nkpt,EIGR2D.nband),dtype=complex)
  fan_add = N.array(zeros((len(temp_info),EIGR2D.nkpt,EIGR2D.nband),dtype=complex))
  ddw_add = N.array(zeros((len(temp_info),EIGR2D.nkpt,EIGR2D.nband),dtype=complex))

  bose = get_bose(EIGR2D.natom,omega,temp_info)

# Get reduced displacement (scaled with frequency)
  displ_red_FAN2,displ_red_DDW2 = get_reduced_displ(EIGR2D.natom,eigvect,omega,gprimd)
# Einstein sum make the vector matrix multiplication ont the correct indices
  fan_corrQ = N.einsum('ijklmn,olnkm->oij',EIGR2D.EIG2D,displ_red_FAN2)
  ddw_corrQ = N.einsum('ijklmn,olnkm->oij',ddw_save,displ_red_DDW2)

  fan_corr = N.einsum('ijk,il->ljk',fan_corrQ,2*bose+1.0)
  ddw_corr = N.einsum('ijk,il->ljk',ddw_corrQ,2*bose+1.0)

  print "Now compute active space ..."

# Now compute active space
  temp = N.einsum('ijklm,nlokp->ijmnpo',GKK,displ_red_FAN2)
  fan_addQ = N.einsum('ijklmn,ijmnk->ijkl',temp,N.conjugate(GKK))
  temp = N.einsum('ijklm,nlokp->ijmnpo',ddw_save2,displ_red_DDW2)
  ddw_addQ = N.einsum('ijklmn,ijmnk->ijkl',temp,N.conjugate(ddw_save2))

  if type == 1: 
    occtmp = EIGR2D.occ[0,0,:]/2 # jband
    delta_E_ddw = N.einsum('ij,k->ijk',eig0[0,:,:].real,N.ones(EIGR2D.nband)) - \
              N.einsum('ij,k->ikj',eig0[0,:,:].real,N.ones(EIGR2D.nband)) - \
              N.einsum('ij,k->ijk',N.ones((EIGR2D.nkpt,EIGR2D.nband)),(2*occtmp-1))*smearing*1j

    tmp = N.einsum('ijkl,lm->mijk',ddw_addQ,2*bose+1.0) # tmp,ikpt,iband,jband
    ddw_add = N.einsum('ijkl,jkl->ijk',tmp,1.0/delta_E_ddw)
    delta_E = N.einsum('ij,k->ijk',eig0[0,:,:].real,N.ones(EIGR2D.nband)) - \
              N.einsum('ij,k->ikj',eigq.EIG[0,:,:].real,N.ones(EIGR2D.nband)) - \
              N.einsum('ij,k->ijk',N.ones((EIGR2D.nkpt,EIGR2D.nband)),(2*occtmp-1))*smearing*1j # ikpt,iband,jband
    omegatmp = omega[:].real # imode
    num1 = N.einsum('ij,k->ijk',bose,N.ones(EIGR2D.nband)) +1.0 \
          - N.einsum('ij,k->ijk',N.ones((3*EIGR2D.natom,len(temp_info))),occtmp) #imode,tmp,jband
    deno1 = N.einsum('ijk,l->ijkl',delta_E,N.ones(3*EIGR2D.natom)) \
          - N.einsum('ijk,l->ijkl',N.ones((EIGR2D.nkpt,EIGR2D.nband,EIGR2D.nband)),omegatmp) #ikpt,iband,jband,imode
    div1 = N.einsum('ijk,lmki->ijklm',num1,1.0/deno1) # (imode,tmp,jband)/(ikpt,iband,jband,imode) ==> imode,tmp,jband,ikpt,iband
    num2 = N.einsum('ij,k->ijk',bose,N.ones(EIGR2D.nband)) \
          + N.einsum('ij,k->ijk',N.ones((3*EIGR2D.natom,len(temp_info))),occtmp) #imode,tmp,jband
    deno2 = N.einsum('ijk,l->ijkl',delta_E,N.ones(3*EIGR2D.natom)) \
          + N.einsum('ijk,l->ijkl',N.ones((EIGR2D.nkpt,EIGR2D.nband,EIGR2D.nband)),omegatmp) #ikpt,iband,jband,imode
    div2 = N.einsum('ijk,lmki->ijklm',num2,1.0/deno2) # (imode,tmp,jband)/(ikpt,iband,jband,imode) ==> imode,tmp,jband,ikpt,iband
    fan_add = N.einsum('ijkl,lmkij->mij',fan_addQ,div1+div2) # ikpt,iband,jband,imode

    print "Now compute generalized g2F Eliashberg electron-phonon spectral function ..."

    fan_tmp = N.sum(fan_addQ,axis=2) # (ikpt, iband, imode)
    ddw_tmp = N.sum(ddw_addQ,axis=2) # (ikpt, iband, imode)
    g_kk = fan_tmp - ddw_tmp

    # Generalized Eliashberg function with Gaussian broadening
    a2F =  zeros((len(energy),EIGR2D.nkpt,EIGR2D.nband),dtype=complex)
    PDOS = zeros((len(energy)),dtype=complex)
    ii = 0
    for ifreq in energy:
      for imode in N.arange(3*EIGR2D.natom):
        a2F[ii,:,:] += g_kk[:,:,imode]*gaussian(ifreq,omegatmp[imode],gaussian_smearing)
        PDOS[ii] += gaussian(ifreq,omegatmp[imode],gaussian_smearing)
      ii += 1


# Correction from active space 
  fan_corr += fan_add
  ddw_corr += ddw_add
  eigen_corr = (fan_corr[:,:,:] - ddw_corr[:,:,:])*wtq
# Index 0 contains total ZPR
  total_corr[0,:,:,:] = eigen_corr[:,:,:]
# Index 1 contains the FAN correction
  total_corr[1,:,:,:] = fan_corr[:,:,:]*wtq
# Index 2 contains the DDW correction
  total_corr[2,:,:,:] = ddw_corr[:,:,:]*wtq
# Index 3-energy contains the Eliashberg functions
  total_corr[3:len(energy)+3,0,:,:] = a2F*wtq

  total_corr = make_average(EIGR2D.nkpt,EIGR2D.nband,degen,total_corr,temp=True)

# Index energy+4 contains the PDOS
  total_corr[len(energy)+4:len(energy)*2+4,0,0,0] = PDOS*wtq

  return total_corr

#########################################################################################################
# Compute total weigth
def compute_wtq(arguments,type):
  if type == 1:
    nbqpt,wtq,eigq_files,DDB_files,EIGR2D_files,GKK_files = arguments
  EIGR2D = system(directory='.',filename=EIGR2D_files)
# If the calculation is on a Homogenous q-point mesh
# retreve the weight of the q-point
  if (wtq == 0):
    wtq = EIGR2D.wtq
    wtq = wtq[0]
  return wtq


class zpm:
  total_corr = None
  def __init__(self,arguments,ddw_save,ddw_save2,nb_cpus,type,temperature,temp_info,smearing,eig0,degen,energy,gaussian_smearing):
    # Parallelize the work over cpus
    pool = multiprocessing.Pool(processes=nb_cpus)
    partial_compute_wtq = partial(compute_wtq,type=type)
    total = pool.imap(partial_compute_wtq,arguments)
    self.total_wtq = sum(total)

# TYPE 1
    if (type == 1 and not temperature):
      partial_dynamic_zpm =  partial(dynamic_zpm,ddw_save=ddw_save,ddw_save2=ddw_save2,type=type,smearing=smearing,eig0=eig0,\
                             degen=degen,energy=energy,gaussian_smearing=gaussian_smearing)
      total = pool.imap(partial_dynamic_zpm,arguments)
      self.total_corr = sum(total)

    if (type == 1 and temperature):
      partial_dynamic_zpm_temp =  partial(dynamic_zpm_temp,ddw_save=ddw_save,ddw_save2=ddw_save2,type=type,temp_info=temp_info,\
                                   smearing=smearing,eig0=eig0,degen=degen,energy=energy,gaussian_smearing=gaussian_smearing)
      total = pool.imap(partial_dynamic_zpm_temp,arguments)
      self.total_corr = sum(total)

    pool.close()
    pool.join()

