
import numpy as np
from numpy import zeros
import netCDF4 as nc


def convert_ddb_to_netcdf(DDB_fname, DDB_nc_fname):
    """Read a _DDB file and convert it into a _DDB.nc file."""
    converter = DdbFileConverter()
    converter.read_txt(DDB_fname)
    converter.write_nc(DDB_nc_fname)


class DdbFileConverter(object):
    """Class to read a DDB file in txt format and write it in netCDF format."""
    # Note that in the netCDF file, it is written as
    # E2D[iat,icart,jat,jcart]
    # While it is read from the txt file as
    # E2D[icart,iat,jcart,jat]

    def read_nc(self, fname):
        """Open the DDB.nc file and read it."""

        with nc.Dataset(fname, 'r') as root:

            self.natom = len(root.dimensions['number_of_atoms'])
            self.ncart = len(root.dimensions['number_of_cartesian_directions'])  # 3
            self.ntypat = len(root.dimensions['number_of_atom_species'])
            #self.nkpt = len(root.dimensions['number_of_kpoints'])  # Not relevant
            #self.nband = len(root.dimensions['max_number_of_states'])  # Not even there
            #self.nsppol = len(root.dimensions['number_of_spins'])  # Not relevant

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

    def write_nc(self, fname):
        """Open withe a DDB.nc file."""

        with nc.Dataset(fname, 'w') as root:

            root.createDimension('number_of_atoms', self.natom)
            root.createDimension('number_of_cartesian_directions', 3)
            root.createDimension('number_of_atom_species', self.ntypat)
            root.createDimension('cplex', 2)

            data = root.createVariable('atom_species', 'i', ('number_of_atoms'))
            data[...] = self.typat[...]

            data = root.createVariable('atomic_masses_amu', 'd', ('number_of_atom_species'))
            data[...] = self.amu[...]

            data = root.createVariable('primitive_vectors', 'd', ('number_of_cartesian_directions', 'number_of_cartesian_directions'))
            data[...] = self.rprim[...]

            data = root.createVariable('reduced_atom_positions', 'd', ('number_of_atoms', 'number_of_cartesian_directions'))
            data[...] = self.xred[...]

            data = root.createVariable('q_point_reduced_coord', 'd', ('number_of_cartesian_directions'))
            data[...] = self.qred[...]

            data = root.createVariable('born_effective_charge_tensor', 'd', ('number_of_cartesian_directions',
                                                                            'number_of_atoms', 'number_of_cartesian_directions'))
            data[...] = self.BECT[...]

            # E2D dimensions must be swapped
            E2D_swap = np.zeros((self.natom, self.ncart, self.natom, self.ncart), dtype=np.complex)
            for iat in range(self.natom):
                for icart in range(self.ncart):
                    for jat in range(self.natom):
                        for jcart in range(self.ncart):
                            E2D_swap[iat, icart, jat, jcart] = self.E2D[icart, iat, jcart, jat]

            data = root.createVariable('second_derivative_of_energy', 'd', ('number_of_atoms', 'number_of_cartesian_directions',
                                                                            'number_of_atoms', 'number_of_cartesian_directions',
                                                                            'cplex'))
            E2D_swap = np.einsum('aibj->bjai', E2D_swap)  # Indicies are reversed when writing them from Fortran.
            data[...,0] = E2D_swap.real
            data[...,1] = E2D_swap.imag


    def read_txt(self, fname):
      """Open the DDB file and read it."""
      self.nkpt=0
      self.rprim = zeros((3,3))
      with open(fname,'r') as DDB:
        Flag = 0
        Flag2 = False
        Flag3 = False
        ikpt = 0
        for line in DDB:
          if line.find('natom') > -1:
            self.natom = np.int(line.split()[1])
          if line.find('nkpt') > -1:
            self.nkpt = np.int(line.split()[1])
            self.kpt  = zeros((self.nkpt,3))
          if line.find('ntypat') > -1:
            self.ntypat = np.int(line.split()[1])
          if line.find('nband') > -1:
            self.nband = np.int(line.split()[1])
          if line.find('acell') > -1:
            line = line.replace('D','E')
            tmp = line.split()
            self.acell = [np.float(tmp[1]),np.float(tmp[2]),np.float(tmp[3])]
          if Flag2:
            line = line.replace('D','E')
            for ii in np.arange(3,self.ntypat):
              self.amu[ii] = np.float(line.split()[ii-3])
              Flag2 = False
          if line.find('amu') > -1:
            line = line.replace('D','E')
            self.amu = zeros((self.ntypat))
            if self.ntypat > 3:
              for ii in np.arange(3):
                self.amu[ii] = np.float(line.split()[ii+1])
                Flag2 = True 
            else:
              for ii in np.arange(self.ntypat):
                self.amu[ii] = np.float(line.split()[ii+1])
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
            for ii in np.arange(12,self.natom): 
              self.typat[ii] = np.float(line.split()[ii-12]) 
            Flag3 = False 
          if line.find(' typat') > -1:
            self.typat = zeros((self.natom))
            if self.natom > 12:
              for ii in np.arange(12):
                self.typat[ii] = np.float(line.split()[ii+1])
                Flag3 = True
            else:
              for ii in np.arange(self.natom):
                self.typat[ii] = np.float(line.split()[ii+1])
          # Read the actual d2E/dRdR matrix
          if Flag == 3:
            line = line.replace('D','E')
            tmp = line.split()
            if not tmp:
              break
            self.E2D[int(tmp[0])-1,int(tmp[1])-1,int(tmp[2])-1,int(tmp[3])-1] = \
              complex(float(tmp[4]),float(tmp[5]))
          # Read the current Q-point
          if line.find('qpt') > -1:
            line = line.replace('D','E')
            tmp = line.split()
            #self.qred = [np.float(tmp[1]),np.float(tmp[2]),np.float(tmp[3])]
            self.qred = np.array([np.float(tmp[1]),np.float(tmp[2]),np.float(tmp[3])])
            Flag = 3
            self.E2D = zeros((3,self.natom,3,self.natom),dtype=complex)

        self.ncart = 3
        self.BECT = zeros((self.ncart,self.natom,self.ncart), dtype=float)
        self.xred = zeros((self.natom, self.ncart), dtype=float)
        for i in range(3):
            self.rprim[i,:] *= self.acell[i]
        self.acell = np.ones(3, dtype=float)
