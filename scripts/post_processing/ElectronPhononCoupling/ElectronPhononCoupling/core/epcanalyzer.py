from __future__ import print_function

import warnings
from copy import copy

import numpy as np
import netCDF4 as nc

from .constants import Ha2eV

from .util import create_directory, formatted_array_lines

from .qptanalyzer import QptAnalyzer

from .mpi import MPI, comm, size, rank, master_only, mpi_watch, i_am_master

# =========================================================================== #

__all__ = ['EpcAnalyzer']


class EpcAnalyzer(object):
    """
    Main class for analysing electron-phonon coupling related quantities.

    It is intented to analyse the files produced by ABINIT
    in a phonon response-function calculation, with one q-point per dataset,
    the first q-point being Gamma.

    """
    verbose = False
    wtq = None
    broadening = None
    temperatures = []
    omegase = []
    smearing = None
    mu = None

    qred = None
    omega = None

    zero_point_renormalization = None
    zero_point_broadening = None
    temperature_dependent_renormalization = None
    temperature_dependent_broadening = None
    zero_point_renormalization_modes = None


    renormalization_is_dynamical = False
    broadening_is_dynamical = False

    self_energy = None
    spectral_function = None
    self_energy_T = None
    spectral_function_T = None

    my_iqpts = [0]

    def __init__(self,
                 nqpt=1,
                 wtq=[1.0],
                 eig0_fname='',
                 eigq_fnames=list(),
                 DDB_fnames=list(),
                 EIGR2D_fnames=list(),
                 EIGI2D_fnames=list(),
                 FAN_fnames=list(),
                 GKK_fnames=list(),
                 output='epc.out',
                 temp_range=[0,0,1],
                 omega_range=[0,0,1],
                 smearing=0.00367,
                 asr=True,
                 fermi_level = None,
                 verbose=False,
                 **kwargs):

        # Check that the minimum number of files is present
        if not eig0_fname:
            raise Exception('Must provide a file for eig0_fname')

        if not EIGR2D_fnames:
            raise Exception('Must provide at least one file for EIGR2D_fnames')

        if not DDB_fnames:
            raise Exception('Must provide at least one file for DDB_fnames')

        if len(wtq) != nqpt:
            raise Exception("Must provide nqpt weights in the 'wtq' list.")

        # Set basic quantities
        self.nqpt = nqpt
        self.set_weights(wtq)

        # Set file names
        self.eig0_fname = eig0_fname
        self.eigq_fnames = eigq_fnames
        self.DDB_fnames = DDB_fnames
        self.EIGR2D_fnames = EIGR2D_fnames
        self.EIGI2D_fnames = EIGI2D_fnames
        self.FAN_fnames = FAN_fnames
        self.GKK_fnames = GKK_fnames

        # Initialize a single QptAnalyzer
        self.qptanalyzer = QptAnalyzer(
            wtq=self.wtq[0],
            eig0_fname=self.eig0_fname,
            DDB_fname=self.DDB_fnames[0],
            EIGR2D0_fname=self.EIGR2D_fnames[0],
            FAN0_fname=self.FAN_fnames[0] if self.FAN_fnames else None,
            GKK0_fname=self.GKK_fnames[0] if self.GKK_fnames else None,
            asr=asr,
            )

        # Read the first DDB and check that it is Gamma
        self.check_gamma()

        # Read other files at q=0 and broadcast the data
        self.read_zero_files()

        # Split the workload between workers
        self.distribute_workload()

        # Get arrays dimensions
        self.nkpt = self.qptanalyzer.eigr2d0.nkpt
        self.nband = self.qptanalyzer.eigr2d0.nband
        self.natom = self.qptanalyzer.eigr2d0.natom
        self.kpts = self.qptanalyzer.eigr2d0.kpt[:,:]

        # Set parameters
        self.set_temp_range(temp_range)
        self.set_omega_range(omega_range)
        self.set_smearing(smearing)
        self.set_output(output)
        self.set_fermi_level(fermi_level)

        self.verbose = verbose

    @master_only
    def check_gamma(self):
        self.qptanalyzer.read_nonzero_files()
        if not self.qptanalyzer.is_gamma:
            raise Exception('The first Q-point is not Gamma.')

    @mpi_watch
    def read_zero_files(self):
        """Read the q=0 files and broadcast to all mpi workers."""

        # Master reads the files
        if i_am_master:
            self.qptanalyzer.read_zero_files()

        # Broadcast
        self.qptanalyzer.broadcast_zero_files()

    def set_temp_range(self, temp_range=(0, 0, 1)):
        """Set the minimum, makimum and step temperature."""
        self.temperatures = np.arange(*temp_range, dtype=float)
        self.ntemp = len(self.temperatures)
        self.qptanalyzer.temperatures = self.temperatures

    def check_temperatures(self):
        if not len(self.temperatures):
            warnings.warn('Temperatures were not set. '
                          'Please specify it with the "temp_range" '
                          'keyword argument ')

    def set_omega_range(self, omega_range=(0, 0, 1)):
        """Set the minimum, makimum and step frequency for the self-energy."""
        self.omegase = np.arange(*omega_range, dtype=float)
        self.nomegase = len(self.omegase)
        self.qptanalyzer.omegase = self.omegase

    def set_smearing(self, smearing_Ha):
        """Set the smearing, in Hartree."""
        self.smearing = smearing_Ha
        self.qptanalyzer.smearing = smearing_Ha
    
    def set_output(self, root):
        """Set the root for output names."""
        self.output = root

    def set_weights(self, wtq, normalize=True):
        """Set the q-points weights."""
        if normalize:
            self.wtq = np.array(wtq) / sum(wtq)
        else:
            self.wtq = np.array(wtq)

    def set_fermi_level(self, mu):
        """Set the Fermi level, in Hartree."""
        self.mu = mu
        self.qptanalyzer.mu = mu

    @mpi_watch
    def find_fermi_level(self):
        """
        Find the Fermi level by gathering information from all workers
        and broadcast the result.
        """
        all_max_val = self.gather_qpt_function('get_max_val')
        all_min_cond = self.gather_qpt_function('get_min_cond')
        if i_am_master:
            max_val = np.max(all_max_val)
            min_cond = np.min(all_min_cond)
            mu = (max_val + min_cond) / 2.0
            mu = np.array(mu, dtype=np.float64)
        else:
            mu = np.empty(1, dtype=np.float64)

        comm.Bcast([mu, MPI.DOUBLE])

        self.set_fermi_level(mu)

    def set_iqpt(self, iqpt):
        """
        Give the qptanalyzer the weight and files corresponding
        to one particular qpoint and read the files. 
        """
        self.qptanalyzer.wtq = self.wtq[iqpt]
        self.qptanalyzer.ddb.fname = self.DDB_fnames[iqpt]

        if self.EIGR2D_fnames:
            self.qptanalyzer.eigr2d.fname = self.EIGR2D_fnames[iqpt]

        if self.eigq_fnames:
            self.qptanalyzer.eigq.fname = self.eigq_fnames[iqpt]

        if self.FAN_fnames:
            self.qptanalyzer.fan.fname = self.FAN_fnames[iqpt]

        if self.GKK_fnames:
            self.qptanalyzer.gkk.fname = self.GKK_fnames[iqpt]

        if self.EIGI2D_fnames:
            self.qptanalyzer.eigi2d.fname = self.EIGI2D_fnames[iqpt]

        self.qptanalyzer.read_nonzero_files()

    def set_ddb(self, iqpt):
        """
        Give the qptanalyzer the weight and ddb file corresponding
        to one particular qpoint, then read and diagonalize the dynamical matrix.
        """
        self.qptanalyzer.wtq = self.wtq[iqpt]
        self.qptanalyzer.ddb.fname = self.DDB_fnames[iqpt]
        self.qptanalyzer.read_ddb()

    @mpi_watch
    def distribute_workload(self):
        """Distribute the q-points indicies to be treated by each worker."""

        max_nqpt_per_worker = self.nqpt // size + min(self.nqpt % size, 1)
        n_active_workers = self.nqpt // max_nqpt_per_worker + min(self.nqpt % max_nqpt_per_worker, 1)

        self.my_iqpts = list()

        for i in range(max_nqpt_per_worker):

            iqpt = rank * max_nqpt_per_worker + i

            if iqpt >= self.nqpt:
                break

            self.my_iqpts.append(iqpt)

    @property
    def active_worker(self):
        return bool(self.my_iqpts)

    def get_active_ranks(self):
        """Get the ranks of all active workers."""
        max_nqpt_per_worker = self.nqpt // size + min(self.nqpt % size, 1)
        n_active_workers = self.nqpt // max_nqpt_per_worker + min(self.nqpt % max_nqpt_per_worker, 1)
        return np.arange(n_active_workers)

    @mpi_watch
    def sum_qpt_function(self, func_name, *args, **kwargs):
        """Call a certain function or each q-points and sum the result."""

        partial_sum = self.sum_qpt_function_me(func_name, *args, **kwargs)

        if i_am_master:
            total = partial_sum
            active_ranks = self.get_active_ranks()
            if len(active_ranks) > 1:
                for irank in active_ranks[1:]:
                    partial_sum = comm.recv(source=irank, tag=irank)
                    total += partial_sum

        elif self.active_worker:
            comm.send(partial_sum, dest=0, tag=rank)
            return

        else:
            return

        # Now I could broadcast the total result to all workers
        # but right now there is no need to.

        return total

    def sum_qpt_function_me(self, func_name, *args, **kwargs):
        """Call a certain function or each q-points of this worker and sum the result."""
        if not self.active_worker:
            return None

        iqpt = self.my_iqpts[0]
        self.set_iqpt(iqpt)

        if self.verbose:
            print("Q-point: {} with wtq = {} and reduced coord. {}".format(
                  iqpt, self.qptanalyzer.wtq, self.qptanalyzer.qred))

        q0 = getattr(self.qptanalyzer, func_name)(*args, **kwargs)
        total = copy(q0)

        if len(self.my_iqpts) == 1:
            return total

        for iqpt in self.my_iqpts[1:]:

            self.set_iqpt(iqpt)

            if self.verbose:
                print("Q-point: {} with wtq = {} and reduced coord. {}".format(
                      iqpt, self.qptanalyzer.wtq, self.qptanalyzer.qred))


            qpt = getattr(self.qptanalyzer, func_name)(*args, **kwargs)
            total += qpt

        return total

    @mpi_watch
    def gather_qpt_function(self, func_name, *args, **kwargs):
        """Call a certain function or each q-points and gather all results."""

        partial = self.gather_qpt_function_me(func_name, *args, **kwargs)

        if i_am_master:
            total = np.zeros([self.nqpt] + list(partial.shape[1:]), dtype=partial.dtype)
            for i, arr in enumerate(partial):
                total[i,...] = arr[...]

            active_ranks = self.get_active_ranks()
            if len(active_ranks) > 1:
                for irank in active_ranks[1:]:
                    partial = comm.recv(source=irank, tag=irank)
                    for arr in partial:
                        i += 1
                        total[i,...] = arr[...]

        elif self.active_worker:
            comm.send(partial, dest=0, tag=rank)
            return

        else:
            return

        # Now I could broadcast the total result to all workers
        # but right now there is no need to.

        return total

    def gather_qpt_function_me(self, func_name, *args, **kwargs):
        """Call a certain function or each q-points of this worker and gather all results."""
        if not self.active_worker:
            return None

        nqpt_me = len(self.my_iqpts)

        iqpt = self.my_iqpts[0]
        self.set_iqpt(iqpt)

        if self.verbose:
            print("Q-point: {} with wtq = {} and reduced coord. {}".format(
                  iqpt, self.qptanalyzer.wtq, self.qptanalyzer.qred))

        q0 = np.array(getattr(self.qptanalyzer, func_name)(*args, **kwargs))
        total = np.zeros([nqpt_me] + list(q0.shape), dtype=q0.dtype)

        if len(self.my_iqpts) == 1:
            return total

        for i, iqpt in enumerate(self.my_iqpts[1:]):

            self.set_iqpt(iqpt)

            if self.verbose:
                print("Q-point: {} with wtq = {} and reduced coord. {}".format(
                      iqpt, self.qptanalyzer.wtq, self.qptanalyzer.qred))


            qpt = getattr(self.qptanalyzer, func_name)(*args, **kwargs)
            total[i+1,...] = qpt[...]

        return total

    @mpi_watch
    def gather_qpt_info(self):
        """Gather qpt reduced coordinates and mode frequencies."""

        partial = self.gather_qpt_info_me()

        if i_am_master:

            qred_all = np.zeros((self.nqpt, 3), dtype=np.float)
            omega_all = np.zeros((self.nqpt, 3 * self.natom), dtype=np.float)

            qred_p, omega_p = partial
            for i, (qred, omega) in enumerate(zip(qred_p, omega_p)):
                qred_all[i,...] = qred[...]
                omega_all[i,...] = omega[...]

            active_ranks = self.get_active_ranks()
            if len(active_ranks) > 1:
                for irank in active_ranks[1:]:
                    partial = comm.recv(source=irank, tag=10000+irank)
                    qred_p, omega_p = partial
                    for qred, omega in zip(qred_p, omega_p):
                        i += 1
                        qred_all[i,...] = qred[...]
                        omega_all[i,...] = omega[...]

        elif self.active_worker:
            comm.send(partial, dest=0, tag=10000+rank)
            return
        else:
            return

        self.qred = qred_all
        self.omega = omega_all

        return self.qred, self.omega

    def gather_qpt_info_me(self):
        """Gather qpt reduced coordinates and mode frequencies."""
        if not self.active_worker:
            return None

        nqpt_me = len(self.my_iqpts)

        qred = np.zeros((nqpt_me, 3), dtype=np.float)
        omega = np.zeros((nqpt_me, 3 * self.natom), dtype=np.float)

        for i, iqpt in enumerate(self.my_iqpts):

            self.set_ddb(iqpt)
            qred[i,:] = self.qptanalyzer.qred[:]
            omega[i,:] = np.real(self.qptanalyzer.omega[:])

        return qred, omega

    def compute_static_zp_renormalization(self):
        """Compute the zero-point renormalization in a static scheme."""
        self.zero_point_renormalization = self.sum_qpt_function('get_zpr_static')
        self.renormalization_is_dynamical = False

    def compute_static_td_renormalization(self):
        """
        Compute the temperature-dependent renormalization in a static scheme.
        """
        self.check_temperatures()
        self.temperature_dependent_renormalization = self.sum_qpt_function('get_tdr_static')
        self.renormalization_is_dynamical = False

    def compute_dynamical_td_renormalization(self):
        """
        Compute the temperature-dependent renormalization in a dynamical scheme.
        """
        self.check_temperatures()
        self.temperature_dependent_renormalization = self.sum_qpt_function('get_tdr_dynamical')
        self.renormalization_is_dynamical = True

    def compute_dynamical_zp_renormalization(self):
        """Compute the zero-point renormalization in a dynamical scheme."""
        self.zero_point_renormalization = self.sum_qpt_function('get_zpr_dynamical')
        self.renormalization_is_dynamical = True

    def compute_static_control_td_renormalization(self):
        """
        Compute the temperature-dependent renormalization in a static scheme
        with the transitions split between active and sternheimer.
        """
        self.check_temperatures()
        self.temperature_dependent_renormalization = self.sum_qpt_function('get_tdr_static_active')
        self.renormalization_is_dynamical = False

    def compute_static_control_zp_renormalization(self):
        """
        Compute the zero-point renormalization in a static scheme
        with the transitions split between active and sternheimer.
        """
        self.zero_point_renormalization = self.sum_qpt_function('get_zpr_static_active')
        self.renormalization_is_dynamical = False

    def compute_static_td_broadening(self):
        """
        Compute the temperature-dependent broadening in a static scheme
        from the EIGI2D files.
        """
        self.check_temperatures()
        self.temperature_dependent_broadening = self.sum_qpt_function('get_tdb_static')
        self.broadening_is_dynamical = False

    def compute_static_zp_broadening(self):
        """
        Compute the zero-point broadening in a static scheme
        from the EIGI2D files.
        """
        self.zero_point_broadening = self.sum_qpt_function('get_zpb_static')
        self.broadening_is_dynamical = False

    def compute_dynamical_td_broadening(self):
        self.check_temperatures()
        warnings.warn("Dynamical lifetime at finite temperature is not yet implemented...proceed with static lifetime")
        return self.compute_static_td_broadening()

    def compute_dynamical_zp_broadening(self):
        """
        Compute the zero-point broadening in a dynamical scheme.
        """
        self.zero_point_broadening = self.sum_qpt_function('get_zpb_dynamical')
        self.broadening_is_dynamical = True

    def compute_static_control_td_broadening(self):
        self.check_temperatures()
        warnings.warn("Static lifetime at finite temperature with control over smearing is not yet implemented...proceed with static lifetime")
        return self.compute_static_td_broadening()

    def compute_static_control_zp_broadening(self):
        """
        Compute the zero-point broadening in a static scheme
        from the FAN files.
        """
        self.zero_point_broadening = self.sum_qpt_function('get_zpb_static_active')
        self.broadening_is_dynamical = False

    def compute_static_control_zp_renormalization_modes(self):
        """
        Compute the zero-point renormalization in a static scheme
        with the transitions split between active and sternheimer.
        Retain the mode decomposition of the zpr.
        """

        self.zero_point_renormalization_modes = self.sum_qpt_function('get_zpr_static_active_modes')
        self.renormalization_is_dynamical = False

    def compute_zp_self_energy(self):
        """
        Compute the zp frequency-dependent self-energy from one q-point.
    
        The self-energy is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is
    
            Simga'_kn(omega) = Sigma_kn(omega + E^0_kn)
    
        """
        self.self_energy = self.sum_qpt_function('get_zp_self_energy')

    def compute_td_self_energy(self):
        """
        Compute the td frequency-dependent self-energy from one q-point.
    
        The self-energy is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is
    
            Simga'_kn(omega,T) = Sigma_kn(omega + E^0_kn,T)
    
        """
        # Make sure the fermi level is set
        if self.mu is None:
            self.find_fermi_level()

        self.self_energy_T = self.sum_qpt_function('get_td_self_energy')

    @master_only
    def compute_zp_spectral_function(self):
        """
        Compute the spectral function of all quasiparticles in the
        semi-static approximation, that is, the 'upper bands' contribution
        to the self-energy is evaluated at the bare energy.

        The spectral function is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is

            A'_kn(omega) = A_kn(omega + E^0_kn)

        """
        nomegase = self.nomegase
        nkpt = self.nkpt
        nband = self.nband

        self.spectral_function = np.zeros((nomegase, nkpt, nband), dtype=float)

        omega = np.einsum('kn,l->knl', np.ones((nkpt, nband)), self.omegase)

        self.spectral_function = (
            (1 / np.pi) * np.abs(self.self_energy.imag) /
            ((omega - self.self_energy.real) ** 2 + self.self_energy.imag ** 2)
            )

    @master_only
    def compute_td_spectral_function(self):
        """
        Compute the spectral function of all quasiparticles in the
        semi-static approximation, that is, the 'upper bands' contribution
        to the self-energy is evaluated at the bare energy.

        The spectral function is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is

            A'_kn(omega) = A_kn(omega + E^0_kn)

        """
        nomegase = self.nomegase
        nkpt = self.nkpt
        nband = self.nband
        ntemp = self.ntemp

        self.spectral_function_T = np.zeros((nomegase, ntemp, nkpt, nband),
                                            dtype=float)

        omega = np.einsum('ijt,l->ijlt', np.ones((nkpt, nband, ntemp)), self.omegase)

        self.spectral_function_T = (
            (1 / np.pi) * np.abs(self.self_energy_T.imag) /
            ((omega - self.self_energy_T.real) ** 2 + self.self_energy_T.imag ** 2)
            )


    @master_only
    def write_netcdf(self):
        """Write all data to a netCDF file."""
        fname = str(self.output) + '_EP.nc'
        create_directory(fname)

        # Write on a NC files with etsf-io name convention
        ncfile = nc.Dataset(fname, 'w')

        # Read dim from first EIGR2D file
        root = nc.Dataset(self.EIGR2D_fnames[0], 'r')

        # Determine nsppol from reading occ
        nsppol = len(root.variables['occupations'][:,0,0])
        if nsppol > 1:
          warnings.warn("nsppol > 1 has not been tested.")
        mband = len(root.dimensions['product_mband_nsppol']) / nsppol

        # Create dimension
        ncfile.createDimension('number_of_atoms',len(root.dimensions['number_of_atoms']))
        ncfile.createDimension('number_of_kpoints',len(root.dimensions['number_of_kpoints']))
        ncfile.createDimension('product_mband_nsppol',len(root.dimensions['product_mband_nsppol']))
        ncfile.createDimension('cartesian',3)
        ncfile.createDimension('cplex',2)
        ncfile.createDimension('number_of_qpoints', self.nqpt)
        ncfile.createDimension('number_of_spins',len(root.dimensions['number_of_spins']))
        ncfile.createDimension('max_number_of_states',mband)
        ncfile.createDimension('number_of_modes', 3 * len(root.dimensions['number_of_atoms']))

        ncfile.createDimension('number_of_temperature',len(self.temperatures))
        ncfile.createDimension('number_of_frequencies',len(self.omegase))

        # Create variable
        data = ncfile.createVariable('reduced_coordinates_of_kpoints','d',
                                     ('number_of_kpoints','cartesian'))
        data[:,:] = root.variables['reduced_coordinates_of_kpoints'][:,:]
        data = ncfile.createVariable('eigenvalues','d',
                                     ('number_of_spins','number_of_kpoints','max_number_of_states'))
        data[:,:,:] = root.variables['eigenvalues'][:,:,:]
        data = ncfile.createVariable('occupations','i',
                                     ('number_of_spins','number_of_kpoints','max_number_of_states'))
        data[:,:,:] = root.variables['occupations'][:,:,:]
        data = ncfile.createVariable('primitive_vectors','d',('cartesian','cartesian'))
        data[:,:] = root.variables['primitive_vectors'][:,:]

        root.close()

        data = ncfile.createVariable('renormalization_is_dynamical', 'i1')
        data[:] = self.renormalization_is_dynamical

        data = ncfile.createVariable('broadening_is_dynamical', 'i1')
        data[:] = self.broadening_is_dynamical

        data = ncfile.createVariable('temperatures','d',('number_of_temperature'))
        data[:] = self.temperatures[:]

        data = ncfile.createVariable('smearing', 'd')
        data[:] = self.smearing

        data = ncfile.createVariable('omegase','d',('number_of_frequencies'))
        data[:] = self.omegase[:]

        zpr = ncfile.createVariable('zero_point_renormalization','d',
            ('number_of_spins', 'number_of_kpoints', 'max_number_of_states'))

        #fan = ncfile.createVariable('fan_zero_point_renormalization','d',
        #    ('number_of_spins', 'number_of_kpoints', 'max_number_of_states'))

        #ddw = ncfile.createVariable('ddw_zero_point_renormalization','d',
        #    ('number_of_spins', 'number_of_kpoints', 'max_number_of_states'))

        if self.zero_point_renormalization is not None:
            zpr[0,:,:] = self.zero_point_renormalization[:,:].real  # FIXME number of spin
            #fan[0,:,:] = self.fan_zero_point_renormalization[:,:].real
            #ddw[0,:,:] = self.ddw_zero_point_renormalization[:,:].real

        data = ncfile.createVariable('temperature_dependent_renormalization','d',
            ('number_of_spins','number_of_kpoints', 'max_number_of_states','number_of_temperature'))

        if self.temperature_dependent_renormalization is not None:
            data[0,:,:,:] = self.temperature_dependent_renormalization[:,:,:].real  # FIXME number of spin

        data = ncfile.createVariable('zero_point_broadening','d',
            ('number_of_spins', 'number_of_kpoints', 'max_number_of_states'))

        if self.zero_point_broadening is not None:
            data[0,:,:] = self.zero_point_broadening[:,:].real  # FIXME number of spin

        data = ncfile.createVariable('temperature_dependent_broadening','d',
            ('number_of_spins','number_of_kpoints', 'max_number_of_states','number_of_temperature'))

        if self.temperature_dependent_broadening is not None:
            data[0,:,:,:] = self.temperature_dependent_broadening[:,:,:].real  # FIXME number of spin

        self_energy = ncfile.createVariable('self_energy','d',
            ('number_of_spins', 'number_of_kpoints', 'max_number_of_states',
             'number_of_frequencies', 'cplex'))

        if self.self_energy is not None:
            self_energy[0,:,:,:,0] = self.self_energy[:,:,:].real  # FIXME number of spin
            self_energy[0,:,:,:,1] = self.self_energy[:,:,:].imag  # FIXME number of spin

        self_energy_T = ncfile.createVariable('self_energy_temperature_dependent','d',
            ('number_of_spins', 'number_of_kpoints', 'max_number_of_states',
             'number_of_frequencies', 'number_of_temperature', 'cplex'))

        if self.self_energy_T is not None:
            self_energy_T[0,:,:,:,:,0] = self.self_energy_T[:,:,:,:].real  # FIXME number of spin
            self_energy_T[0,:,:,:,:,1] = self.self_energy_T[:,:,:,:].imag  # FIXME number of spin

        spectral_function = ncfile.createVariable('spectral_function','d',
            ('number_of_spins', 'number_of_kpoints', 'max_number_of_states',
             'number_of_frequencies'))

        if self.spectral_function is not None:
            spectral_function[0,:,:,:] = self.spectral_function[:,:,:]  # FIXME number of spin

        spectral_function_T = ncfile.createVariable('spectral_function_temperature_dependent','d',
            ('number_of_spins', 'number_of_kpoints', 'max_number_of_states',
             'number_of_frequencies', 'number_of_temperature'))

        if self.spectral_function_T is not None:
            spectral_function_T[0,:,:,:,:] = self.spectral_function_T[:,:,:,:]  # FIXME number of spin

        zpr_modes = ncfile.createVariable('zero_point_renormalization_by_modes','d',
            ('number_of_modes', 'number_of_spins', 'number_of_kpoints', 'max_number_of_states'))
        if self.zero_point_renormalization_modes is not None:
            zpr_modes[:,0,:,:] = self.zero_point_renormalization_modes[:,:,:]

        data = ncfile.createVariable('reduced_coordinates_of_qpoints','d',
                                     ('number_of_qpoints', 'cartesian'))
        if self.qred is not None:
            data[...] = self.qred[...]

        data = ncfile.createVariable('phonon_mode_frequencies','d',
                                     ('number_of_qpoints', 'number_of_modes'))
        if self.omega is not None:
            data[...] = self.omega[...]

        ncfile.close()


    @master_only
    def write_renormalization(self):
        """Write the computed renormalization in a text file."""
        fname = str(self.output) + "_REN.txt"
        create_directory(fname)

        with open(fname, "w") as O:

            if self.zero_point_renormalization is not None:
                O.write("Total zero point renormalization (eV) for {} Q points\n".format(self.nqpt))
                for ikpt, kpt in enumerate(self.kpts):
                    O.write('Kpt: {0[0]} {0[1]} {0[2]}\n'.format(kpt))
                    for line in formatted_array_lines(self.zero_point_renormalization[ikpt,:].real*Ha2eV):
                        O.write(line)

            if self.temperature_dependent_renormalization is not None:
                O.write("Temperature dependence at Gamma (eV)\n")
                for iband in range(self.nband):
                  O.write('Band: {}\n'.format(iband))
                  for tt, T in enumerate(self.temperatures):
                    ren = self.temperature_dependent_renormalization[0,iband,tt].real * Ha2eV
                    O.write("{:>8.1f}  {:>12.8f}\n".format(T, ren))


    @master_only
    def write_broadening(self):
        """Write the computed broadening in a text file."""
        fname = str(self.output) + "_BRD.txt"
        create_directory(fname)

        with open(fname, "w") as O:

            if self.zero_point_broadening is not None:
                O.write("Total zero point broadening (eV) for {} Q points\n".format(self.nqpt))
                for ikpt, kpt in enumerate(self.kpts):
                    O.write('Kpt: {0[0]} {0[1]} {0[2]}\n'.format(kpt))
                    for line in formatted_array_lines(self.zero_point_broadening[ikpt,:].real*Ha2eV):
                        O.write(line)

            if self.temperature_dependent_broadening is not None:
                O.write("Temperature dependence at Gamma\n")
                for iband in range(self.nband):
                  O.write('Band: {}\n'.format(iband))
                  for tt, T in enumerate(self.temperatures):
                    brd = self.temperature_dependent_broadening[0,iband,tt].real * Ha2eV
                    O.write("{:>8.1f}  {:>12.8f}\n".format(T, brd))

