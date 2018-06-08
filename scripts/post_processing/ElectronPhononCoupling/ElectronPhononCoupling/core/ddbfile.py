
from __future__ import print_function

__author__ = "Gabriel Antonius, Samuel Ponce"

import os
import warnings

import numpy as np
from numpy import zeros
import netCDF4 as nc

from .mpi import MPI, comm, size, rank, mpi_watch

from .constants import tol5, tol6, me_amu, kb_HaK
from .functions import get_bose
from . import EpcFile

__all__ = ['DdbFile']


class DdbFile(EpcFile):

    _rprim = np.identity(3)
    gprimd = np.identity(3)
    omega = None
    asr = True
    
    def __init__(self, *args, **kwargs):
        self.asr = kwargs.pop('asr', True)
        super(DdbFile, self).__init__(*args, **kwargs)

    def set_amu(self, amu):
        """
        Set the values for the atom masses.

        Arguments:
            amu: [ntypat]
                Atom masses for each atom type, in atomic mass units.
        """
        self.amu = np.array(amu)

    def read_nc(self, fname=None):
        """Open the DDB.nc file and read it."""
        fname = fname if fname else self.fname

        super(DdbFile, self).read_nc(fname)

        with nc.Dataset(fname, 'r') as root:

            self.natom = len(root.dimensions['number_of_atoms'])
            self.ncart = len(root.dimensions['number_of_cartesian_directions'])  # 3
            self.ntypat = len(root.dimensions['number_of_atom_species'])

            self.typat = root.variables['atom_species'][:self.natom]
            self.amu = root.variables['atomic_masses_amu'][:self.ntypat]
            self.rprim = root.variables['primitive_vectors'][:self.ncart,:self.ncart]
            self.xred = root.variables['reduced_atom_positions'][:self.natom,:self.ncart]
            self.qred = root.variables['q_point_reduced_coord'][:]

            # The d2E/dRdR' matrix
            self.E2D = np.zeros((self.natom, self.ncart, self.natom, self.ncart), dtype=np.complex)
            self.E2D.real = root.variables['second_derivative_of_energy'][:,:,:,:,0]
            self.E2D.imag = root.variables['second_derivative_of_energy'][:,:,:,:,1]
            self.E2D = np.einsum('aibj->bjai', self.E2D)  # Indicies are reversed when writing them from Fortran.

            self.BECT = root.variables['born_effective_charge_tensor'][:self.ncart,:self.natom,:self.ncart]

    def broadcast(self):
        """Broadcast the data from master to all workers."""
    
        comm.Barrier()

        if rank == 0:
            dim = np.array([self.natom, self.ncart, self.ntypat], dtype=np.int)
        else:
            dim = np.empty(3, dtype=np.int)

        comm.Bcast([dim, MPI.INT])

        if rank != 0:

            self.natom, self.ncart, self.ntypat = dim[:]

            self.typat = np.empty(self.natom, dtype='i')
            self.amu = np.empty(self.ntypat, dtype=np.float)
            rprim = np.empty((self.ncart, self.ncart), dtype=np.float)
            self.xred = np.empty((self.natom, self.ncart), dtype=np.float)
            self.qred = np.empty((self.ncart), dtype=np.float)

            self.E2D = np.empty((self.natom, self.ncart, self.natom, self.ncart),
                                dtype=np.complex)

            self.BECT = np.empty((self.ncart, self.natom, self.ncart), dtype=np.float)

        else:
            rprim = self.rprim

        comm.Bcast([self.typat, MPI.INT])
        comm.Bcast([self.amu, MPI.DOUBLE])
        comm.Bcast([self.xred, MPI.DOUBLE])
        comm.Bcast([self.qred, MPI.DOUBLE])
        comm.Bcast([self.E2D, MPI.COMPLEX])
        comm.Bcast([self.BECT, MPI.DOUBLE])
        comm.Bcast([rprim, MPI.DOUBLE])

        self.rprim = rprim

    @property
    def is_gamma(self):
        return np.allclose(self.qred,[0.0,0.0,0.0])

    @property
    def rprim(self):
        return self._rprim

    @rprim.setter
    def rprim(self, value):
        self._rprim = np.array(value)
        self.gprimd = np.linalg.inv(np.matrix(self._rprim))

    @property
    def nmode(self):
        return 3 * self.natom

    def get_mass_scaled_dynmat_cart(self):
        """
        Format the dynamical matrix in a 3Nx3N matrix,
        scale with masses, and transform into Cartesian coordinates.
        """
        # Retrive the amu for each atom
        amu = zeros(self.natom)
        for ii in np.arange(self.natom):
          jj = self.typat[ii]
          amu[ii] = self.amu[jj-1]
    
        # Transform from 2nd-order matrix (non-cartesian coordinates, 
        # masses not included, asr not included ) from self to
        # dynamical matrix, in cartesian coordinates, asr not imposed.
        E2D_cart = zeros((3,self.natom,3,self.natom),dtype=complex)
        for ii in np.arange(self.natom):
          for jj in np.arange(self.natom):
            for dir1 in np.arange(3):
              for dir2 in np.arange(3):
                for dir3 in np.arange(3):
                  for dir4 in np.arange(3):
                    E2D_cart[dir1,ii,dir2,jj] += (self.E2D[ii,dir3,jj,dir4] *
                            self.gprimd[dir1,dir3] * self.gprimd[dir2,dir4])

        # Reduce the 4 dimensional E2D_cart matrice to 2 dimensional Dynamical matrice
        # with scaled masses.
        Dyn_mat = zeros((3*self.natom,3*self.natom),dtype=complex)
        for ii in np.arange(self.natom):
          for dir1 in np.arange(3):
            ipert1 = ii * 3 + dir1
            for jj in np.arange(self.natom):
              for dir2 in np.arange(3):
                ipert2 = jj * 3 + dir2

                Dyn_mat[ipert1,ipert2] = (E2D_cart[dir1,ii,dir2,jj] *
                                          me_amu / np.sqrt(amu[ii]*amu[jj]))
    
        # Hermitianize the dynamical matrix
        #dynmat = np.matrix(Dyn_mat)
        dynmat = Dyn_mat
        dynmat = 0.5 * (dynmat + dynmat.transpose().conjugate())

        return dynmat


    def get_E2D_cart(self):
        """
        Transform the 2nd-order matrix from non-cartesian coordinates
        to cartesian coordinates...and also swap atom and cartesian indicies.
        """
        # Transform from 2nd-order matrix (non-cartesian coordinates, 
        # masses not included, asr not included ) from self to
        # dynamical matrix, in cartesian coordinates, asr not imposed.
        E2D_cart = zeros((3,self.natom,3,self.natom),dtype=complex)
        for ii in np.arange(self.natom):
          for jj in np.arange(self.natom):
            for dir1 in np.arange(3):
              for dir2 in np.arange(3):
                for dir3 in np.arange(3):
                  for dir4 in np.arange(3):
                    E2D_cart[dir1,ii,dir2,jj] += (self.E2D[ii,dir3,jj,dir4] *
                            self.gprimd[dir1,dir3] * self.gprimd[dir2,dir4])

        return E2D_cart
                                        

    def compute_dynmat(self, asr=None, zero_negative=True):
        """
        Diagonalize the dynamical matrix.
    
        Returns:
          omega: the frequencies, in Ha
          eigvect: the eigenvectors, in reduced coord
        """
        asr = asr if asr is not None else self.asr
    
        # Retrive the amu for each atom
        amu = zeros(self.natom)
        for ii in np.arange(self.natom):
          jj = self.typat[ii]
          amu[ii] = self.amu[jj-1]
    
        dynmat = self.get_mass_scaled_dynmat_cart()
    
        # Diagonalize the matrix
        eigval, eigvect = np.linalg.eigh(dynmat)
    
        # Scale the eigenvectors 
        for ii in np.arange(self.natom):
          for dir1 in np.arange(3):
            ipert = ii * 3 + dir1
            eigvect[ipert] = eigvect[ipert] * np.sqrt(me_amu / amu[ii])

        # Nullify imaginary frequencies
        if zero_negative:
            for i, eig in enumerate(eigval):
              if eig < 0.0:
                warnings.warn("An eigenvalue is negative with value: {} ... but proceed with value 0.0".format(jj))
                eigval[i] = 0.0

        # Impose the accoustic sum rule
        if asr and self.is_gamma:
          eigval[0] = 0.0
          eigval[1] = 0.0
          eigval[2] = 0.0

        # Frequencies
        self.omega = np.sqrt(np.abs(eigval)) * np.sign(eigval)
        self.eigvect = eigvect
    
        return self.omega, self.eigvect

    def get_reduced_displ(self):
        """
        Compute the mode eigenvectors, scaled by the mode displacements
        Also transform from cartesian to reduced coordinates.

        Returns: polvec[nmode,3,natom]
        """

        # Minimal value for omega (Ha)
        omega_tolerance = 1e-5

        self.polvec = zeros((self.nmode,3,self.natom), dtype=complex)
        xi_at = zeros(3, dtype=complex)

        omega, eigvect = self.compute_dynmat()

        for imode in range(self.nmode):
            # Skip mode with zero frequency (leave displacements null)
            if omega[imode].real < omega_tolerance:
              continue

            z0 = 1. / np.sqrt(2.0 * omega[imode].real)

            for iatom in np.arange(self.natom):
                for idir in range(3):
                    xi_at[idir] = eigvect[3*iatom+idir,imode] * z0

                for idir in range(3):
                    for jdir in range(3):
                        self.polvec[imode,idir,iatom] += xi_at[jdir] * self.gprimd[jdir,idir]

        return self.polvec
    
    def get_reduced_displ_squared(self):
        """
        Compute the squared reduced displacements (scaled by phonon frequencies)
        for the Fan and the DDW terms.
        """
        # Minimal value for omega (Ha)
        omega_tolerance = 1e-5

        natom = self.natom
        omega, eigvect = self.compute_dynmat()

        displ_FAN = zeros((3,3), dtype=complex)
        displ_DDW = zeros((3,3), dtype=complex)
        displ_red_FAN2 = zeros((3*natom,natom,natom,3,3), dtype=complex)
        displ_red_DDW2 = zeros((3*natom,natom,natom,3,3), dtype=complex)

        for imode in np.arange(3*natom):

          # Skip mode with zero frequency (leave displacements null)
          if omega[imode].real < omega_tolerance:
            continue

          for iatom1 in np.arange(natom):
            for iatom2 in np.arange(natom):
              for idir1 in np.arange(0,3):
                for idir2 in np.arange(0,3):

                  displ_FAN[idir1,idir2] = (
                      eigvect[3*iatom2+idir2,imode].conj() *
                      eigvect[3*iatom1+idir1,imode] / (2.0 * omega[imode].real)
                      )

                  displ_DDW[idir1,idir2] = (
                      eigvect[3*iatom2+idir2,imode].conj() *
                      eigvect[3*iatom2+idir1,imode] +
                      eigvect[3*iatom1+idir2,imode].conj() *
                      eigvect[3*iatom1+idir1,imode]) / (4.0 * omega[imode].real)

              # Now switch to reduced coordinates in 2 steps (more efficient)
              tmp_displ_FAN = zeros((3,3),dtype=complex)
              tmp_displ_DDW = zeros((3,3),dtype=complex)

              for idir1 in np.arange(3):
                for idir2 in np.arange(3):

                  tmp_displ_FAN[:,idir1] = tmp_displ_FAN[:,idir1] + displ_FAN[:,idir2] * self.gprimd[idir2,idir1]
                  tmp_displ_DDW[:,idir1] = tmp_displ_DDW[:,idir1] + displ_DDW[:,idir2] * self.gprimd[idir2,idir1]

              displ_red_FAN = zeros((3,3),dtype=complex)
              displ_red_DDW = zeros((3,3),dtype=complex)
              for idir1 in np.arange(3):
                for idir2 in np.arange(3):
                  displ_red_FAN[idir1,:] = displ_red_FAN[idir1,:] + tmp_displ_FAN[idir2,:] * self.gprimd[idir2,idir1]
                  displ_red_DDW[idir1,:] = displ_red_DDW[idir1,:] + tmp_displ_DDW[idir2,:] * self.gprimd[idir2,idir1]
    
              displ_red_FAN2[imode,iatom1,iatom2,:,:] = displ_red_FAN[:,:]
              displ_red_DDW2[imode,iatom1,iatom2,:,:] = displ_red_DDW[:,:]

        self.displ_red_FAN2 = displ_red_FAN2
        self.displ_red_DDW2 = displ_red_DDW2
    
        return displ_red_FAN2, displ_red_DDW2

    def get_bose(self, temperatures):
        """
        Get the Bose-Einstein occupations on a range of temperatures.
        Returns: bose(3*natom, Ntemperatures)
        """
        if self.omega is None:
            self.compute_dynmat()

        bose = zeros((3*self.natom, len(temperatures)))

        #bose[:,:] = get_bose(self.omega, temperatures)
        for imode, omega in enumerate(self.omega):
            bose[imode,:] = get_bose(omega, temperatures)

        return bose


    # This old function reads the DDB from the ascii file.
    # It is left here for legacy.
    #
    #def DDB_file_open(self, filefullpath):
    #  """Open the DDB file and read it."""
    #  if not (os.path.isfile(filefullpath)):
    #    raise Exception('The file "%s" does not exists!' %filefullpath)
    #  with open(filefullpath,'r') as DDB:
    #    Flag = 0
    #    Flag2 = False
    #    Flag3 = False
    #    ikpt = 0
    #    for line in DDB:
    #      if line.find('natom') > -1:
    #        self.natom = np.int(line.split()[1])
    #      if line.find('nkpt') > -1:
    #        self.nkpt = np.int(line.split()[1])
    #        self.kpt  = zeros((self.nkpt,3))
    #      if line.find('ntypat') > -1:
    #        self.ntypat = np.int(line.split()[1])
    #      if line.find('nband') > -1:
    #        self.nband = np.int(line.split()[1])
    #      if line.find('acell') > -1:
    #        line = line.replace('D','E')
    #        tmp = line.split()
    #        self.acell = [np.float(tmp[1]),np.float(tmp[2]),np.float(tmp[3])]
    #      if Flag2:
    #        line = line.replace('D','E')
    #        for ii in np.arange(3,self.ntypat):
    #          self.amu[ii] = np.float(line.split()[ii-3])
    #          Flag2 = False
    #      if line.find('amu') > -1:
    #        line = line.replace('D','E')
    #        self.amu = zeros((self.ntypat))
    #        if self.ntypat > 3:
    #          for ii in np.arange(3):
    #            self.amu[ii] = np.float(line.split()[ii+1])
    #            Flag2 = True 
    #        else:
    #          for ii in np.arange(self.ntypat):
    #            self.amu[ii] = np.float(line.split()[ii+1])
    #      if line.find(' kpt ') > -1:
    #        line = line.replace('D','E')
    #        tmp = line.split()
    #        self.kpt[0,0:3] = [float(tmp[1]),float(tmp[2]),float(tmp[3])]
    #        ikpt = 1
    #        continue
    #      if ikpt < self.nkpt and ikpt > 0:
    #        line = line.replace('D','E')
    #        tmp = line.split()
    #        self.kpt[ikpt,0:3] = [float(tmp[0]),float(tmp[1]),float(tmp[2])]  
    #        ikpt += 1
    #        continue
    #      if Flag == 2:
    #        line = line.replace('D','E')
    #        tmp = line.split()
    #        self.rprim[2,0:3] = [float(tmp[0]),float(tmp[1]),float(tmp[2])]
    #        Flag = 0
    #      if Flag == 1:
    #        line = line.replace('D','E')
    #        tmp = line.split()
    #        self.rprim[1,0:3] = [float(tmp[0]),float(tmp[1]),float(tmp[2])]
    #        Flag = 2
    #      if line.find('rprim') > -1:
    #        line = line.replace('D','E')
    #        tmp = line.split()
    #        self.rprim[0,0:3] = [float(tmp[1]),float(tmp[2]),float(tmp[3])]
    #        Flag = 1
    #      if Flag3:
    #        line = line.replace('D','E')
    #        for ii in np.arange(12,self.natom): 
    #          self.typat[ii] = np.float(line.split()[ii-12]) 
    #        Flag3 = False 
    #      if line.find(' typat') > -1:
    #        self.typat = zeros((self.natom))
    #        if self.natom > 12:
    #          for ii in np.arange(12):
    #            self.typat[ii] = np.float(line.split()[ii+1])
    #            Flag3 = True
    #        else:
    #          for ii in np.arange(self.natom):
    #            self.typat[ii] = np.float(line.split()[ii+1])
    #      # Read the actual d2E/dRdR matrix
    #      if Flag == 3:
    #        line = line.replace('D','E')
    #        tmp = line.split()
    #        if not tmp:
    #          break
    #        self.E2D[int(tmp[0])-1,int(tmp[1])-1,int(tmp[2])-1,int(tmp[3])-1] = \
    #          complex(float(tmp[4]),float(tmp[5]))
    #      # Read the current Q-point
    #      if line.find('qpt') > -1:
    #        line = line.replace('D','E')
    #        tmp = line.split()
    #        self.iqpt = [np.float(tmp[1]),np.float(tmp[2]),np.float(tmp[3])]
    #        Flag = 3
    #        self.E2D = zeros((3,self.natom,3,self.natom),dtype=complex)


