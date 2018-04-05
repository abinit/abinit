from __future__ import print_function

__author__ = "Gabriel Antonius"

from copy import copy

import numpy as np
import netCDF4 as nc

from . import EpcFile

from .constants import tol6, kb_HaK

from .mpi import MPI, comm, size, rank, mpi_watch

__all__ = ['EigFile']


class EigFile(EpcFile):

    def __init__(self, *args, **kwargs):
        super(EigFile, self).__init__(*args, **kwargs)
        self.EIG = None
        self.degen = None

    def read_nc(self, fname=None):
        """Open the Eig.nc file and read it."""
        fname = fname if fname else self.fname

        super(EigFile, self).read_nc(fname)

        with nc.Dataset(fname, 'r') as root:

            # nspin, nkpt, nband
            self.EIG = root.variables['Eigenvalues'][:,:,:] 

            # nkpt, 3
            self.Kptns = root.variables['Kptns'][:,:]

    @mpi_watch
    def broadcast(self):
        """Broadcast the data from master to all workers."""
        comm.Barrier()

        if rank == 0:
            dim = np.array([self.nspin, self.nkpt, self.nband], dtype=np.int)
        else:
            dim = np.empty(3, dtype=np.int)

        comm.Bcast([dim, MPI.INT])

        if rank != 0:
            self.EIG = np.empty(dim, dtype=np.float64)
            self.Kptns = np.empty((dim[1],3), dtype=np.float64)

        comm.Bcast([self.EIG, MPI.DOUBLE])
        comm.Bcast([self.Kptns, MPI.DOUBLE])

    @property
    def nspin(self):
        return self.EIG.shape[0] if self.EIG is not None else None

    @property
    def nkpt(self):
        return self.EIG.shape[1] if self.EIG is not None else None

    @property
    def nband(self):
        return self.EIG.shape[2] if self.EIG is not None else None

    def iter_spin_band_eig(self, ikpt):
        """
        Iterator over spin index, band index, and eigenvalues at one k-point.
        Yields tuples (ispin, iband, eig) in order of increasing eig.
        """
        nspin, nkpt, nband = self.EIG.shape
        sbe = [(ispin, 0, self.EIG[ispin,ikpt,0]) for ispin in range(nspin)]
        cmp_sbe = lambda sbe1, sbe2: cmp(sbe1[2], sbe2[2])
        while sbe:
            min_sbe = sorted(sbe, cmp=cmp_sbe)[0]
            yield min_sbe
    
            i = sbe.index(min_sbe)
            ispin, iband, eig = min_sbe
    
            if iband == nband - 1:
                del sbe[i]
            else:
                sbe[i] = (ispin, iband+1, self.EIG[ispin, ikpt, iband+1])
    
    def get_degen(self):
        """
        Compute the degeneracy of the bands.
    
        Returns
        -------
        degen: 2D list (nkpt, )
            For each k-point, contains a list of groups of (s, n) tuples
            which are degenerated.
            For example, if there is only a single spin (spin unpolarized case)
            and two k-points, and at the first k-point there are two triplets,
            then the output is
                [[[(0,1), (0,2), (0,3)],  [(0,4), (0,5), (0,6)]], []]
    
        """
        nspin, nkpt, nband = self.EIG.shape
    
        degen = list()
        for ikpt in range(nkpt):
    
            kpt_degen = list()
            group = list()
            last_ispin, last_iband, last_eig = 0, 0, -float('inf')
    
            for sbe in self.iter_spin_band_eig(ikpt):
                ispin, iband, eig = sbe
    
                if np.isclose(last_eig, eig, rtol=1e-12, atol=1e-5):
                    if not group:
                        group.append((last_ispin, last_iband))
                    group.append((ispin, iband))
    
                else:
                    if group:
                        kpt_degen.append(group)
                        group = list()
    
                last_ispin, last_iband, last_eig = ispin, iband, eig
    
            degen.append(kpt_degen)

        self.degen = degen
    
        return degen

    def make_average(self, arr):
        """ 
        Average a quantity over degenerated states.
        Does not work with spin yet.
    
        Arguments
        ---------
    
        arr: numpy.ndarray(..., nkpt, nband)
            An array of any dimension, of which the two last indicies are
            the kpoint and the band.
    
        Returns
        -------
    
        arr: numpy.ndarray(..., nkpt, nband)
            The array with the values of degenerated bands averaged.
    
        """

        if not self.degen:
            self.get_degen()

        nkpt, nband = arr.shape[-2:]
    
        for ikpt in range(nkpt):
            for group in self.degen[ikpt]:
                average = copy(arr[...,ikpt,group[0][1]])
                for ispin, iband in group[1:]:
                    average += arr[...,ikpt,iband]
    
                average /= len(group)
                for ispin, iband in group:
                    arr[...,ikpt,iband] = average
    
        return arr

    def symmetrize_fan_degen(self, fan_epc):
        """
        Enforce coupling terms to be zero on the diagonal
        and in degenerate states subset.
    
        Arguments
        ---------
        fan_epc: np.ndarray, shape=(nkpt,nband,nband,nmode)
            The coupling matrix V_ij V_ji
    
        Returns
        -------
    
        fan_epc_symmetrized: np.ndarray, shape=(nkpt,nband,nband,nmode)
        """
        if not self.degen:
            self.get_degen()

        nkpt, nband, mband, nmode = fan_epc.shape
     
        offdiag = np.zeros((nkpt, nband, nband))
        offdiag[:] =  np.ones((nband, nband)) - np.identity(nband)
    
        for ikpt in range(nkpt):
            for group in self.degen[ikpt]:
                for degi in group:
                    for degj in group:
                        ieig, jeig = degi[1], degj[1]
                        offdiag[ikpt][ieig][jeig] = 0
    
        fan_epc_sym = np.einsum('ijkl,ijk->ijkl', fan_epc, offdiag)
    
        return fan_epc_sym

    def get_fermi_function_T0(self, mu):
        """
        Get the Fermi function for T=0.
        Returns: occ[nspin,nkpt,nband]
        """
        occ = np.zeros((self.nspin, self.nkpt, self.nband))
        occ[np.where(self.EIG < mu)] = 1.0
        occ[np.where(self.EIG > mu)] = 0.0
        return occ

    def get_fermi_function(self, mu, temperatures):
        """
        Compute the Fermi function for the occupations,
        given a chemical potential.
        """

        ntemp = len(temperatures)
        occ = np.zeros((self.nspin, self.nkpt, self.nband, ntemp))

        for itemp, T in enumerate(temperatures):

            if T < tol6:
                occ[...,itemp] = self.get_fermi_function_T0(mu)
                continue

            beta = 1. / (kb_HaK * T)

            #occ[...,itemp] = 1. / (np.exp(beta * (self.EIG - mu)) + 1)

            for ispin in range(self.nspin):
                for ikpt in range(self.nkpt):
                    for iband in range(self.nband):

                        betaE = beta * (self.EIG[ispin,ikpt,iband] - mu)

                        if betaE < -20:
                            occ[ispin,ikpt,iband,itemp] = 1
                        elif betaE > 20:
                            occ[ispin,ikpt,iband,itemp] = 0
                        else:
                            occ[ispin,ikpt,iband,itemp] = 1./(np.exp(betaE)+1)

        return occ
 
