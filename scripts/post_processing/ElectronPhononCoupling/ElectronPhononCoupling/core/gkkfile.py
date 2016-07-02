
from __future__ import print_function

__author__ = "Gabriel Antonius"

import os

import numpy as np
from numpy import zeros, einsum
import netCDF4 as nc

from .mpi import MPI, comm, size, rank, mpi_watch

from . import EpcFile

__all__ = ['GkkFile']


class GkkFile(EpcFile):
    
    def read_nc(self, fname=None):
        """Open the EIG2D.nc file and read it."""
        fname = fname if fname else self.fname

        with nc.Dataset(fname, 'r') as root:

            self.natom = len(root.dimensions['number_of_atoms'])
            self.nkpt = len(root.dimensions['number_of_kpoints'])
            self.nband = len(root.dimensions['max_number_of_states'])
            self.nsppol = len(root.dimensions['number_of_spins'])

            # number_of_spins, number_of_kpoints, max_number_of_states
            self.occ = root.variables['occupations'][:,:,:]

            # number_of_spins, number_of_kpoints, max_number_of_states   
            self.eigenvalues = root.variables['eigenvalues'][:,:,:]

            # number_of_kpoints, 3
            self.kpt = root.variables['reduced_coordinates_of_kpoints'][:,:]
            self.qred = root.variables['current_q_point'][:]
            self.wtq = root.variables['current_q_point_weight'][:]
            self.rprimd = root.variables['primitive_vectors'][:,:]

            # nband, natom, ncart, nkpt, product_mband_nsppol*2 
            GKKtmp = root.variables['second_derivative_eigenenergies_actif'][:,:,:,:,:]
            GKKtmp2 = np.einsum('ijkno->nokji', GKKtmp)
            self.GKK = np.zeros((self.nkpt, self.nsppol*self.nband, 3, self.natom, self.nband), dtype=np.complex)
            self.GKK.real[...] = GKKtmp2[:, ::2, ...]
            self.GKK.imag[...] = GKKtmp2[:, 1::2, ...]
            self.GKK = np.reshape(self.GKK,(self.nkpt,self.nsppol,self.nband,3,self.natom,self.nband))

    def get_gkk_squared(self):
        """
        Get squared values of gkk with  reordered indices.
        Returns:
            gkk2[nkpt,nband,3,natom,3,natom,nband]
        """

        # No spin polarization at the moment.

        gkk2 = zeros((self.nkpt, self.nband, 3, self.natom, 3, self.natom,
                      self.nband), dtype=np.complex)

        # nkpt,nband,3,natom,3,natom,nband
        gkk2 = einsum('ijklm,ijnom->ijklnom', self.GKK[:,0,...],
                                              self.GKK[:,0,...].conjugate())

        return gkk2


    @mpi_watch
    def broadcast(self):
        """Broadcast the data from master to all workers."""
    
        comm.Barrier()

        if rank == 0:
            dim = np.array([self.natom, self.nkpt, self.nband, self.nsppol], dtype=np.int)
        else:
            dim = np.empty(4, dtype=np.int)

        comm.Bcast([dim, MPI.INT])

        if rank != 0:

            self.natom, self.nkpt, self.nband, self.nsppol = dim[:]

            self.occ = np.empty((self.nsppol, self.nkpt, self.nband), dtype=np.float)

            self.GKK = np.empty((self.nkpt,self.nsppol,self.nband,3,self.natom,self.nband),
                                dtype=np.complex)

            self.eigenvalues = np.empty((self.nsppol, self.nkpt, self.nband), dtype=np.float)

            # number_of_kpoints, 3
            self.kpt = np.empty((self.nkpt, 3), dtype=np.float)
            self.qred = np.empty(3, dtype=np.float)
            self.wtq = np.empty(1, dtype=np.float)
            self.rprimd = np.empty((3, 3), dtype=np.float)

        comm.Bcast([self.occ, MPI.DOUBLE])
        comm.Bcast([self.GKK, MPI.COMPLEX])
        comm.Bcast([self.eigenvalues, MPI.DOUBLE])
        comm.Bcast([self.kpt, MPI.DOUBLE])
        comm.Bcast([self.qred, MPI.DOUBLE])
        comm.Bcast([self.wtq, MPI.DOUBLE])
        comm.Bcast([self.rprimd, MPI.DOUBLE])
