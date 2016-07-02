#!/bin/bash

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import copy
from matplotlib.colors import LogNorm

class Hexc(object):

    def __init__(self, nsppol, nkpt, nval, ncond, hsize, ham, kpt, vcks2t,lomo, homo,lumo,humo,diag):

        self.nsppol = nsppol
        self.nkpt = nkpt
        self.nval = nval
        self.ncond = ncond
        self.hsize = hsize
        self.ham = ham
        self.vcks2t = vcks2t
        self.kpt = kpt

        self.lomo = lomo
        self.lumo = lumo
        self.homo = homo
        self.humo = humo
	self.diag = diag

    @classmethod
    def from_file(cls,filename):
        root = nc.Dataset(filename, 'r')
        nsppol = len(root.dimensions['number_of_spins'])
        nkpt = len(root.dimensions['number_of_kpoints'])
        nval = len(root.dimensions['max_number_of_valence_bands'])
        ncond = len(root.dimensions['max_number_of_conduction_bands'])
        hsize = len(root.dimensions['total_number_of_transitions'])

        kpt = root.variables['reduced_coordinates_of_kpoints'][:,:] # Nkpt, 3
        ham = root.variables['hamiltonian'][:,:,:] # Itp, It, [0,1]

        ham_c = np.vectorize(complex)(ham[:,:,0],ham[:,:,1])
        vcks2t = root.variables['vcks2t'][:,:,:,:] # S, K, C, V

        lomo = root.variables['lomo'][:]
        homo = root.variables['homo'][:]
        lumo = root.variables['lumo'][:]
        humo = root.variables['humo'][:]
        diag = root.variables['diagonal'][:,:]
        diag_c = np.vectorize(complex)(diag[:,0],diag[:,1])

        return cls(nsppol, nkpt, nval, ncond, hsize, ham_c, kpt, vcks2t, lomo, homo, lumo, humo, diag_c)

    def plot(self,ic=None,iv=None,icp=None,ivp=None,list_ikpt=None, title=None):
        fullmat = None
        chk_bands = (ic is not None and iv is not None and icp is not None and ivp is not None)

        if not chk_bands:
          fullmat = self.ham
        else:
          fullmat = np.zeros((self.nkpt,self.nkpt),dtype=np.complex)
          for is1 in np.arange(self.nsppol):
            for is2 in np.arange(self.nsppol):
              for ik in np.arange(self.nkpt):
                for ikp in np.arange(self.nkpt):
                  it = self.vcks2t[is1,ik,ic-self.lumo[is1],iv-self.lomo[is1]]-1
		  itp = self.vcks2t[is2,ikp,icp-self.lumo[is2],ivp-self.lomo[is2]]-1
                  fullmat[ik,ikp] = self.ham[it,itp] # Fortran -> python

        mymat = None
        if list_ikpt is None:
          mymat = fullmat
        else:
          nkpt2 = len(list_ikpt)
          mymat = np.zeros((nkpt2,nkpt2),dtype=np.complex)
          if chk_bands:
            for ik1 in np.arange(nkpt2):
              for ik2 in np.arange(nkpt2):
                mymat[ik1,ik2] = fullmat[list_ikpt[ik1],list_ikpt[ik2]]
          else:
            mymat = np.zeros((self.nval*self.ncond*nkpt2,self.nval*self.ncond*nkpt2))
            it2 = 0
            for is1 in np.arange(self.nsppol):
              for ik in np.arange(nkpt2):
                for ic in np.arange(self.lumo,self.humo+1):
                  for iv in np.arange(self.lomo,self.homo+1):
                    itp2 = 0
                    for is2 in np.arange(self.nsppol):
                      for ikp in np.arange(nkpt2):
                        for icp in np.arange(self.lumo,self.humo+1):
                          for ivp in np.arange(self.lomo,self.homo+1):
                            it = self.vcks2t[is1,list_ikpt[ik],ic-self.lumo[is1],iv-self.lomo[is1]]-1
                            itp = self.vcks2t[is2,list_ikpt[ikp],icp-self.lumo[is2],ivp-self.lomo[is2]]-1
                            mymat[it2,itp2] = self.ham[it,itp] # Fortran -> python
                            itp2 += 1
                    it2 += 1
        
        fig = plt.figure()
	imgplot = plt.imshow(np.abs(mymat)*self.nkpt,norm=LogNorm(vmin=1e-5,vmax=10),interpolation='none')
        if title is not None:
          plt.title(title)
        plt.colorbar()

    def get_filtered_kpt(self,xval=None,yval=None,zval=None):
       
       list_ikpt = []
    
       for ikpt in np.arange(self.nkpt):
         if (  (xval is None or np.abs(self.kpt[ikpt,0]-xval) < 1e-8) 
           and (yval is None or np.abs(self.kpt[ikpt,1]-yval) < 1e-8)
           and (zval is None or np.abs(self.kpt[ikpt,2]-zval) < 1e-8)):
           list_ikpt.append(ikpt)
    
       return list_ikpt

