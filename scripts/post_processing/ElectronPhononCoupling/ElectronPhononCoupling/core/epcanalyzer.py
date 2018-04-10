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
    Main class for analyzing electron-phonon coupling related quantities.

    It is intented to analyze the files produced by ABINIT
    in a phonon response-function calculation, with one q-point per dataset,
    the first q-point being Gamma.

    For documentation, see `ElectronPhononCoupling.compute`
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

    self_energy = None
    self_energy_T = None
    self_energy_static = None
    self_energy_static_T = None
    self_energy_fan_active = None


    split_fan_ddw = False
    renormalization_is_dynamical = False
    broadening_is_dynamical = False

    self_energy = None
    spectral_function = None
    self_energy_T = None
    spectral_function_T = None

    my_iqpts = [0]

    def __init__(self,

                 # Options
                 asr=True,
                 verbose=False,

                 # Parameters
                 nqpt=1,
                 wtq=[1.0],
                 temp_range=[0,0,1],
                 omega_range=[0,0,1],
                 smearing=0.00367,
                 fermi_level = None,
                 amu = None,

                 double_smearing = False,
                 smearing_width = 0.0367,
                 smearing_above = 0.00367,
                 smearing_below = 0.00367,

                 # File names
                 rootname='epc.out',
                 eigk_fname='',
                 eigq_fnames=list(),
                 ddb_fnames=list(),
                 eigr2d_fnames=list(),
                 eigi2d_fnames=list(),
                 fan_fnames=list(),
                 gkk_fnames=list(),

                 # Double grid
                 nqpt_fine=1,
                 wtq_fine=[1.0],
                 eigq_fine_fnames=list(),
                 gkk_fine_fnames=list(),
                 ddb_fine_fnames=list(),

                 **kwargs):

        # Check that the minimum number of files is present
        if not eigk_fname:
            raise Exception('Must provide a file for eigk_fname')

        if not ddb_fnames:
            raise Exception('Must provide at least one file for ddb_fnames')

        if len(wtq) != nqpt:
            raise Exception("Must provide nqpt weights in the 'wtq' list.")

        # Set basic quantities
        self.nqpt = nqpt
        self.set_weights(wtq)

        self.nqpt_fine = nqpt_fine
        self.set_weights_fine(wtq_fine)

        # Set file names
        self.eig0_fname = eigk_fname
        self.eigq_fnames = eigq_fnames
        self.ddb_fnames = ddb_fnames
        self.eigr2d_fnames = eigr2d_fnames
        self.eigi2d_fnames = eigi2d_fnames
        self.fan_fnames = fan_fnames
        self.gkk_fnames = gkk_fnames

        self.eigq_fine_fnames = eigq_fine_fnames
        self.ddb_fine_fnames = ddb_fine_fnames
        self.gkk_fine_fnames = gkk_fine_fnames


        # Initialize a single QptAnalyzer for q=0

        # Select first gkk file
        if self.gkk_fnames:
            gkk0 = self.gkk_fnames[0]
        elif self.gkk_fine_fnames:
            gkk0 = self.gkk_fine_fnames[0]
        else:
            gkk0 = None

        # Select first fan file
        if self.fan_fnames:
            fan0 = self.fan_fnames[0]
        else:
            fan0 = None

        # Select first eigr2d file
        if self.eigr2d_fnames:
            eigr2d0 = self.eigr2d_fnames[0]
        else:
            eigr2d0 = None

        # Select first ddb file
        if self.ddb_fnames:
            ddb0 = self.ddb_fnames[0]
        else:
            ddb0 = None  # I suppose other things will break...

        self.qptanalyzer = QptAnalyzer(
            wtq=self.wtq[0],
            eigk_fname=self.eig0_fname,
            ddb_fname=ddb0,
            eigr2d0_fname=eigr2d0,
            fan0_fname=fan0,
            gkk0_fname=gkk0,
            asr=asr,
            amu=amu,
            double_smearing = double_smearing,
            smearing_width = smearing_width,
            smearing_above = smearing_above,
            smearing_below = smearing_below,
            )

        # Read the first DDB and check that it is Gamma
        self.check_gamma()

        # Read other files at q=0 and broadcast the data
        self.read_zero_files()

        # Set parameters
        self.set_temp_range(temp_range)
        self.set_omega_range(omega_range)
        self.set_smearing(smearing)
        self.set_rootname(rootname)

        # Split the workload between workers
        # (needed here to find the fermi level)
        self.distribute_workload()

        # Set the fermi level
        if fermi_level is None:
            self.find_fermi_level()
        else:
            self.set_fermi_level(fermi_level)

        self.verbose = verbose

    @property
    def nc_output(self):
        return str(self.rootname) + '_EP.nc'

    @property
    def ren_dat(self):
        return str(self.rootname) + '_REN.dat'

    @property
    def BRD_dat(self):
        return str(self.rootname) + '_BRD.dat'

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

    @property
    def nkpt(self):
        return self.qptanalyzer.nkpt

    @property
    def nband(self):
        return self.qptanalyzer.nband

    @property
    def natom(self):
        return self.qptanalyzer.natom

    @property
    def kpts(self):
        return self.qptanalyzer.kpts

    @property
    def nomegase(self):
        return len(self.omegase)

    @property
    def ntemp(self):
        return len(self.temperatures)

    def set_temp_range(self, temp_range=(0, 0, 1)):
        """Set the minimum, maximum and step temperature."""
        args = list(temp_range)
        assert len(args) == 3
        minimum, maximum, step = args
        if all([isinstance(i, int) for i in args]):
            if (maximum - minimum) % step == 0:
                maximum += 1
        self.temperatures = np.arange(minimum, maximum, step, dtype=float)
        self.qptanalyzer.temperatures = self.temperatures

    def check_temperatures(self):
        if not len(self.temperatures):
            warnings.warn('Temperatures were not set. '
                          'Please specify it with the "temp_range" '
                          'keyword argument ')

    def set_omega_range(self, omega_range=(0, 0, 1)):
        """Set the minimum, makimum and step frequency for the self-energy."""
        self.omegase = np.arange(*omega_range, dtype=float)
        self.qptanalyzer.omegase = self.omegase

    def set_smearing(self, smearing_Ha):
        """Set the smearing, in Hartree."""
        self.smearing = smearing_Ha
        self.qptanalyzer.smearing = smearing_Ha
    
    def set_rootname(self, root):
        """Set the root for output names."""
        self.rootname = root

    def set_weights(self, wtq, normalize=True):
        """Set the q-points weights."""
        if normalize:
            self.wtq = np.array(wtq) / sum(wtq)
        else:
            self.wtq = np.array(wtq)

    def set_weights_fine(self, wtq, normalize=True):
        """Set the q-points weights."""
        if normalize:
            self.wtq_fine = np.array(wtq) / sum(wtq)
        else:
            self.wtq_fine = np.array(wtq)

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

    def set_iqpt(self, iqpt, fine=False):
        """
        Give the qptanalyzer the weight and files corresponding
        to one particular qpoint and read the files. 
        """
        if fine:
            self.set_iqpt_fine(iqpt)
        else:
            self.set_iqpt_coarse(iqpt)

    def set_iqpt_coarse(self, iqpt):
        """
        Give the qptanalyzer the weight and files corresponding
        to one particular qpoint and read the files. 
        """
        self.qptanalyzer.wtq = self.wtq[iqpt]
        self.qptanalyzer.ddb.fname = self.ddb_fnames[iqpt]

        if self.eigr2d_fnames:
            self.qptanalyzer.eigr2d.fname = self.eigr2d_fnames[iqpt]

        if self.eigq_fnames:
            self.qptanalyzer.eigq.fname = self.eigq_fnames[iqpt]

        if self.fan_fnames:
            self.qptanalyzer.fan.fname = self.fan_fnames[iqpt]

        if self.gkk_fnames:
            self.qptanalyzer.gkk.fname = self.gkk_fnames[iqpt]

        if self.eigi2d_fnames:
            self.qptanalyzer.eigi2d.fname = self.eigi2d_fnames[iqpt]

        self.qptanalyzer.read_nonzero_files()

    def set_iqpt_fine(self, iqpt):
        """
        Give the qptanalyzer the weight and files corresponding
        to one particular qpoint and read the files. 
        """
        self.qptanalyzer.wtq = self.wtq_fine[iqpt]
        self.qptanalyzer.ddb.fname = self.ddb_fine_fnames[iqpt]
        self.qptanalyzer.eigq.fname = self.eigq_fine_fnames[iqpt]
        self.qptanalyzer.gkk.fname = self.gkk_fine_fnames[iqpt]
        self.qptanalyzer.read_nonzero_files()

    def set_ddb(self, iqpt):
        """
        Give the qptanalyzer the weight and ddb file corresponding
        to one particular qpoint, then read and diagonalize the dynamical matrix.
        """
        self.qptanalyzer.wtq = self.wtq[iqpt]
        self.qptanalyzer.ddb.fname = self.ddb_fnames[iqpt]
        self.qptanalyzer.read_ddb()

    @mpi_watch
    def distribute_workload(self, fine=False):
        """Distribute the q-points indicies to be treated by each worker."""
        if fine:
            self.my_iqpts = self.get_iqpts(self.nqpt_fine)
        else:
            self.my_iqpts = self.get_iqpts(self.nqpt)

    def get_iqpts(self, nqpt):
        """Distribute the q-points indicies to be treated by each worker."""

        max_nqpt_per_worker = (
            nqpt // size + min(nqpt % size, 1))
        n_active_workers = (
            nqpt // max_nqpt_per_worker
            + min(nqpt % max_nqpt_per_worker, 1))

        my_iqpts = list()

        for i in range(max_nqpt_per_worker):

            iqpt = rank * max_nqpt_per_worker + i

            if iqpt >= nqpt:
                break

            my_iqpts.append(iqpt)

        return my_iqpts

    @property
    def active_worker(self):
        return bool(self.my_iqpts)

    def get_active_ranks(self,fine=False):
        """Get the ranks of all active workers."""
    
        if fine:
           nqpt = self.nqpt_fine
        else:
           nqpt = self.nqpt
 
        #max_nqpt_per_worker = (self.nqpt // size
        #                       + min(self.nqpt % size, 1))
        #n_active_workers = (self.nqpt // max_nqpt_per_worker
        #                    + min(self.nqpt % max_nqpt_per_worker, 1))
        max_nqpt_per_worker = (nqpt // size
                               + min(nqpt % size, 1))
        n_active_workers = (nqpt // max_nqpt_per_worker
                            + min(nqpt % max_nqpt_per_worker, 1))
        return np.arange(n_active_workers)

    @mpi_watch
    def sum_qpt_function(self, func_name, fine=False, *args, **kwargs):
        """Call a certain function or each q-points and sum the result."""

        partial_sum = self.sum_qpt_function_me(func_name, fine=fine,
                                               *args, **kwargs)

        if i_am_master:
            total = partial_sum
            #active_ranks = self.get_active_ranks()
            active_ranks = self.get_active_ranks(fine)
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

    def sum_qpt_function_me(self, func_name, fine=False, *args, **kwargs):
        """
        Call a certain function or each q-points of this worker
        and sum the result.
        """
        if not self.active_worker:
            return None

        iqpt = self.my_iqpts[0]
        self.set_iqpt(iqpt, fine=fine)

        if self.verbose:
            print("Q-point: {} with wtq = {} and reduced coord. {}".format(
                  iqpt, self.qptanalyzer.wtq, self.qptanalyzer.qred))

        q0 = getattr(self.qptanalyzer, func_name)(*args, **kwargs)
        total = copy(q0)

        if len(self.my_iqpts) == 1:
            return total

        for iqpt in self.my_iqpts[1:]:

            self.set_iqpt(iqpt, fine=fine)

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

            # Contruct an array with the shape of partial,
            # adding a dimension of length nqpt.
            total = np.zeros([self.nqpt] + list(partial.shape[1:]),
                             dtype=partial.dtype)

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
        """
        Call a certain function or each q-points of this worker
        and gather all results.
        """
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

        total[0,...] = q0[...]

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

    @mpi_watch
    def sum_qpt_functions_double_grid(self, func_coarse, func_fine,
                                      *args, **kwargs):
        """
        Sum a certain function on the coarse grid,
        and another one on the fine grid.
        Only master sums the result.
        """
        self.distribute_workload(fine=False)
        sum_coarse = self.sum_qpt_function(func_coarse, fine=False)

        self.distribute_workload(fine=True)
        #self.read_zero_files()
        sum_fine = self.sum_qpt_function(func_fine, fine=True)

        if i_am_master:
            result = sum_coarse + sum_fine
        else:
            result = None

        return result

    def compute_static_zp_renormalization_nosplit(self):
        """Compute the zero-point renormalization in a static scheme."""
        self.distribute_workload()
        self.zero_point_renormalization = self.sum_qpt_function('get_zpr_static_sternheimer')
        self.renormalization_is_dynamical = False

    def compute_static_td_renormalization_nosplit(self):
        """
        Compute the temperature-dependent renormalization in a static scheme.
        """
        self.check_temperatures()
        self.distribute_workload()
        self.temperature_dependent_renormalization = self.sum_qpt_function(
            'get_tdr_static_nosplit')
        self.renormalization_is_dynamical = False

    def compute_dynamical_td_renormalization(self):
        """
        Compute the temperature-dependent renormalization
        in a dynamical scheme.
        """
        self.check_temperatures()
        self.distribute_workload()
        self.temperature_dependent_renormalization = self.sum_qpt_function(
            'get_tdr_dynamical')
        self.renormalization_is_dynamical = True

    def compute_dynamical_td_renormalization_double_grid(self):
        """
        Compute the temperature-dependent renormalization
        in a dynamical scheme.
        """
        self.check_temperatures()
        self.temperature_dependent_renormalization = (
            self.sum_qpt_functions_double_grid('get_tdr_static_nosplit',
                                               'get_tdr_dynamical_active'))
        self.renormalization_is_dynamical = True

    def compute_dynamical_zp_renormalization_double_grid(self):
        """
        Compute the temperature-dependent renormalization
        in a dynamical scheme.
        """
        self.check_temperatures()
        self.zero_point_renormalization = (
            self.sum_qpt_functions_double_grid('get_zpr_static_sternheimer',
                                               'get_zpr_dynamical_active'))
        self.renormalization_is_dynamical = True

    def compute_dynamical_zp_renormalization_modes_double_grid(self):
        """
        Compute the temperature-dependent renormalization
        in a dynamical scheme.
        """
        self.check_temperatures()
        self.zero_point_renormalization_modes = (
            self.sum_qpt_functions_double_grid('get_zpr_static_sternheimer_modes',
                                               'get_zpr_dynamical_active_modes'))
        self.renormalization_is_dynamical = True

    def compute_dynamical_zp_renormalization(self):
        """Compute the zero-point renormalization in a dynamical scheme."""
        self.distribute_workload()
        self.zero_point_renormalization = (
            self.sum_qpt_function('get_zpr_dynamical'))
        self.renormalization_is_dynamical = True

    def compute_static_td_renormalization(self):
        """
        Compute the temperature-dependent renormalization in a static scheme
        with the transitions split between active and sternheimer.
        """
        self.check_temperatures()
        self.distribute_workload()
        self.temperature_dependent_renormalization = (
            self.sum_qpt_function('get_tdr_static'))
        self.renormalization_is_dynamical = False

    def compute_static_zp_renormalization(self):
        """
        Compute the zero-point renormalization in a static scheme
        with the transitions split between active and sternheimer.
        """
        self.zero_point_renormalization = (
            self.sum_qpt_function('get_zpr_static'))
        self.renormalization_is_dynamical = False

    def compute_dynamical_td_broadening(self):
        """
        Compute the temperature-dependent broadening in a static scheme
        from the GKK files.
        """
        self.check_temperatures()
        self.distribute_workload()
        self.temperature_dependent_broadening = (
            self.sum_qpt_function('get_tdb_dynamical'))
        self.broadening_is_dynamical = True

    def compute_dynamical_zp_broadening(self):
        """
        Compute the zero-point broadening in a static scheme
        from the GKK files.
        """
        self.distribute_workload()
        self.zero_point_broadening = (
            self.sum_qpt_function('get_zpb_dynamical'))
        self.broadening_is_dynamical = True

    def compute_static_td_broadening(self):
        """
        Compute the temperature-dependent broadening in a static scheme
        from the GKK files.
        """
        self.check_temperatures()
        self.distribute_workload()
        self.temperature_dependent_broadening = (
            self.sum_qpt_function('get_tdb_static'))
        self.broadening_is_dynamical = False

    def compute_static_zp_broadening(self):
        """
        Compute the zero-point broadening in a static scheme
        from the GKK files.
        """
        self.distribute_workload()
        self.zero_point_broadening = (
            self.sum_qpt_function('get_zpb_static'))
        self.broadening_is_dynamical = False

    def compute_static_td_broadening_nosplit(self):
        """
        Compute the temperature-dependent broadening in a static scheme
        from the EIGI2D files.
        """
        self.check_temperatures()
        self.distribute_workload()
        self.temperature_dependent_broadening = (
            self.sum_qpt_function('get_tdb_static_nosplit'))
        self.broadening_is_dynamical = False

    def compute_static_zp_broadening_nosplit(self):
        """
        Compute the zero-point broadening in a static scheme
        from the EIGI2D files.
        """
        self.distribute_workload()
        self.zero_point_broadening = (
            self.sum_qpt_function('get_zpb_static_nosplit'))
        self.broadening_is_dynamical = False

    def compute_ddw_active_zpr(self):
        """
        Compute the zero-point renormalization in a static scheme
        with the transitions split between active and sternheimer.
        """
        self.zero_point_renormalization = (
            self.sum_qpt_function('get_zpr_ddw_active'))
        self.renormalization_is_dynamical = False

    def compute_static_zp_broadening(self):
        """
        Compute the zero-point broadening in a static scheme
        from the FAN files.
        """
        self.distribute_workload()
        self.zero_point_broadening = self.sum_qpt_function('get_zpb_static')
        self.broadening_is_dynamical = False

    def compute_static_zp_renormalization_modes(self):
        """
        Compute the zero-point renormalization in a static scheme
        with the transitions split between active and sternheimer.
        Retain the mode decomposition of the zpr.
        """
        self.distribute_workload()
        self.zero_point_renormalization_modes = (
            self.sum_qpt_function('get_zpr_static_modes'))
        self.renormalization_is_dynamical = False

    def compute_dynamical_zp_renormalization_modes(self):
        """
        Compute the zero-point renormalization in a static scheme
        with the transitions split between active and sternheimer.
        Retain the mode decomposition of the zpr.
        """
        self.distribute_workload()
        self.zero_point_renormalization_modes = (
            self.sum_qpt_function('get_zpr_dynamical_modes'))
        self.renormalization_is_dynamical = True

    def compute_zp_self_energy(self):
        """
        Compute the zero-point frequency-dependent self-energy.
    
        The self-energy is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is
    
            Simga'_kn(omega) = Sigma_kn(omega + E^0_kn)
    
        """
        self.distribute_workload()
        self.self_energy = self.sum_qpt_function('get_zp_self_energy')

    def compute_td_self_energy(self):
        """
        Compute the temperature-dependent frequency-dependent self-energy.
    
        The self-energy is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is
    
            Simga'_kn(omega,T) = Sigma_kn(omega + E^0_kn,T)
    
        """
        self.distribute_workload()
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

        omega = np.einsum('ijt,l->ijlt',
                          np.ones((nkpt, nband, ntemp)), self.omegase)

        self.spectral_function_T = (
            (1 / np.pi) * np.abs(self.self_energy_T.imag) /
            ((omega - self.self_energy_T.real) ** 2
             + self.self_energy_T.imag ** 2)
            )

    def compute_zp_self_energy_double_grid(self):
        """
        Compute the temperature-dependent renormalization
        in a dynamical scheme.
        """
        self.check_temperatures()
        self.self_energy = (
          self.sum_qpt_functions_double_grid('get_zp_self_energy_sternheimer',
                                             'get_zp_self_energy_active'))

    def compute_td_self_energy_double_grid(self):
        """
        Compute the temperature-dependent renormalization
        in a dynamical scheme.
        """
        self.check_temperatures()
        self.self_energy_T = (
          self.sum_qpt_functions_double_grid('get_td_self_energy_sternheimer',
                                             'get_td_self_energy_active'))

    def compute_zp_self_energy_static(self):
        """
        Compute the static part of the zero-point self-energy.
        This includes the Fan and DDW contribution of the Sternheimer space
        and the DDW contribution of the active space.
    
        """
        self.distribute_workload()
        ddw_active = self.sum_qpt_function('get_zpr_ddw_active')
        sternheimer = self.sum_qpt_function('get_zpr_static_sternheimer')

        if i_am_master:
            self.self_energy_static = ddw_active + sternheimer
        else:
            self.self_energy_static = None

        return self.self_energy_static

    def compute_td_self_energy_static(self):
        """
        Compute the static part of the zero-point self-energy.
        This includes the Fan and DDW contribution of the Sternheimer space
        and the DDW contribution of the active space.
        """
        self.distribute_workload()
        ddw_active = self.sum_qpt_function('get_tdr_ddw_active')
        sternheimer = self.sum_qpt_function('get_tdr_static_sternheimer')

        if i_am_master:
            self.self_energy_static_T = ddw_active + sternheimer
        else:
            self.self_energy_static_T = None

        return self.self_energy_static_T
            
    def compute_zp_self_energy_static_double_grid(self):
        """
        Compute the static part of the zero-point self-energy.
        This includes the Fan and DDW contribution of the Sternheimer space
        and the DDW contribution of the active space.
    
        """
        self.self_energy_static = (
            self.sum_qpt_functions_double_grid('get_zpr_static_sternheimer',
                                               'get_zpr_ddw_active'))

        return self.self_energy_static
        
    def compute_td_self_energy_static_double_grid(self):
        """
        Compute the static part of the td self-energy.
        This includes the Fan and DDW contribution of the Sternheimer space
        and the DDW contribution of the active space.
    
        """
        self.self_energy_static_T = (
            self.sum_qpt_functions_double_grid('get_tdr_static_sternheimer',
                                               'get_tdr_ddw_active'))

        return self.self_energy_static_T
    
    def compute_td_self_energy_active(self):
        """
        Compute the temperature-dependent frequency-dependent self-energy
        from the active part only (GKK), neglecting the Sternheimer part
        (EIGR2D).
    
        The self-energy is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is
    
            Simga'_kn(omega,T) = Sigma_kn(omega + E^0_kn,T)
    
        """
        self.distribute_workload()
        self.self_energy_T = self.sum_qpt_function('get_td_self_energy_active')

    def compute_zp_self_energy_active(self):
        """
        Compute the zero-point frequency-dependent self-energy from the active
        part only (GKK), neglecting the Sternheimer part (EIGR2D).
    
        The self-energy is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is
    
            Simga'_kn(omega,T) = Sigma_kn(omega + E^0_kn,T)
    
        """
        self.distribute_workload()
        self.self_energy = self.sum_qpt_function('get_zp_self_energy_active')

    def compute_td_self_energy_static_double_grid(self):
        """
        Compute the static part of the td self-energy.
        This includes the Fan and DDW contribution of the Sternheimer space
        and the DDW contribution of the active space.
    
        """
        self.self_energy_static_T = (
            self.sum_qpt_functions_double_grid('get_tdr_static_sternheimer',
                                               'get_tdr_ddw_active'))

        return self.self_energy_static_T

    def compute_zp_self_energy_fan_active(self):
        """
        Compute only the active part of Fan term for check.

        """
        self.self_energy_fan_active = (
            self.sum_qpt_function('get_zp_fan_active'))

        return self.self_energy_fan_active

    @master_only
    def write_netcdf(self):
        """Write all data to a netCDF file."""

        if self.eigr2d_fnames:
            dim_fname = self.eigr2d_fnames[0]
        elif self.gkk_fnames:
            dim_fname = self.gkk_fnames[0]
        elif self.fan_fnames:
            dim_fname = self.fan_fnames[0]
        else:
            raise Exception('Need at least one file to read the dimensions: ' +
                            'EIGR2D, GKK, or FAN. ' +
                            'How did you even get there?')

        create_directory(self.nc_output)

        # Write on a NC files with etsf-io name convention
        with nc.Dataset(self.nc_output, 'w') as ds:

            # Read dim from first EIGR2D file
            dim = nc.Dataset(dim_fname, 'r')

            # Determine nsppol from reading occ
            nsppol = len(dim.variables['occupations'][:,0,0])
            if nsppol > 1:
              warnings.warn("nsppol > 1 has not been tested.")
            mband = len(dim.dimensions['product_mband_nsppol']) / nsppol

            # Create dimension
            ds.createDimension('number_of_atoms',
                               len(dim.dimensions['number_of_atoms']))
            ds.createDimension('number_of_kpoints',
                               len(dim.dimensions['number_of_kpoints']))
            ds.createDimension('product_mband_nsppol',
                               len(dim.dimensions['product_mband_nsppol']))

            ds.createDimension('cartesian', 3)
            ds.createDimension('cplex', 2)
            ds.createDimension('number_of_qpoints', self.nqpt)
            ds.createDimension('number_of_spins',
                               len(dim.dimensions['number_of_spins']))
            ds.createDimension('max_number_of_states',mband)
            ds.createDimension('number_of_modes',
                               3*len(dim.dimensions['number_of_atoms']))

            ds.createDimension('number_of_temperature', len(self.temperatures))
            ds.createDimension('number_of_frequencies', len(self.omegase))

            # Write data on the eigenvalues
            data = ds.createVariable('reduced_coordinates_of_kpoints', 'd',
                                     ('number_of_kpoints','cartesian'))
            data[:,:] = dim.variables['reduced_coordinates_of_kpoints'][:,:]

            data = ds.createVariable(
                'eigenvalues','d',
                ('number_of_spins','number_of_kpoints','max_number_of_states'))
            data[:,:,:] = dim.variables['eigenvalues'][:,:,:]

            data = ds.createVariable(
                'occupations','i',
                ('number_of_spins','number_of_kpoints','max_number_of_states'))
            data[:,:,:] = dim.variables['occupations'][:,:,:]

            data = ds.createVariable(
                'primitive_vectors', 'd',
                ('cartesian','cartesian'))

            data[:,:] = dim.variables['primitive_vectors'][:,:]

            dim.close()

            # Write epc data
            data = ds.createVariable('renormalization_is_dynamical', 'i1')
            data[:] = self.renormalization_is_dynamical

            data = ds.createVariable('broadening_is_dynamical', 'i1')
            data[:] = self.broadening_is_dynamical

            data = ds.createVariable('temperatures','d',
                                     ('number_of_temperature'))
            data[:] = self.temperatures[:]

            data = ds.createVariable('smearing', 'd')
            data[:] = self.smearing

            data = ds.createVariable('omegase', 'd',
                                     ('number_of_frequencies'))
            data[:] = self.omegase[:]

            # qpt
            data = ds.createVariable(
                'reduced_coordinates_of_qpoints','d',
                ('number_of_qpoints', 'cartesian'))
            if self.qred is not None:
                data[...] = self.qred[...]

            # omega
            data = ds.createVariable(
                'phonon_mode_frequencies','d',
                ('number_of_qpoints', 'number_of_modes'))
            if self.omega is not None:
                data[...] = self.omega[...]


            # ZPR
            zpr = ds.createVariable(
                'zero_point_renormalization','d',
                ('number_of_spins', 'number_of_kpoints',
                 'max_number_of_states'))

            #fan = ds.createVariable(
            #   'fan_zero_point_renormalization','d',
            #   ('number_of_spins', 'number_of_kpoints',
            #    'max_number_of_states'))

            #ddw = ds.createVariable(
            #   'ddw_zero_point_renormalization','d',
            #   ('number_of_spins', 'number_of_kpoints',
            #    'max_number_of_states'))

            if self.zero_point_renormalization is not None:
                # FIXME number of spin
                zpr[0,:,:] = self.zero_point_renormalization[:,:].real
                #fan[0,:,:] = self.fan_zero_point_renormalization[:,:].real
                #ddw[0,:,:] = self.ddw_zero_point_renormalization[:,:].real

            # TDR
            data = ds.createVariable(
                'temperature_dependent_renormalization','d',
                ('number_of_spins','number_of_kpoints',
                 'max_number_of_states','number_of_temperature'))

            if self.temperature_dependent_renormalization is not None:
                # FIXME number of spin
                data[0,:,:,:] = (
                    self.temperature_dependent_renormalization[:,:,:].real)

            # ZPR
            data = ds.createVariable(
                'zero_point_broadening','d',
                ('number_of_spins', 'number_of_kpoints',
                 'max_number_of_states'))

            if self.zero_point_broadening is not None:
                # FIXME number of spin
                data[0,:,:] = self.zero_point_broadening[:,:].real

            zpr_modes = ds.createVariable(
                'zero_point_renormalization_by_modes','d',
                ('number_of_modes', 'number_of_spins', 'number_of_kpoints',
                 'max_number_of_states'))

            if self.zero_point_renormalization_modes is not None:
                zpr_modes[:,0,:,:] = (
                self.zero_point_renormalization_modes[:,:,:])

            # TDB
            data = ds.createVariable(
                'temperature_dependent_broadening','d',
                ('number_of_spins','number_of_kpoints',
                 'max_number_of_states','number_of_temperature'))

            if self.temperature_dependent_broadening is not None:
                # FIXME number of spin
                data[0,:,:,:] = (
                    self.temperature_dependent_broadening[:,:,:].real)

            # ZSE
            self_energy = ds.createVariable('self_energy','d',
                ('number_of_spins', 'number_of_kpoints',
                 'max_number_of_states', 'number_of_frequencies', 'cplex'))

            if self.self_energy is not None:

                # FIXME number of spin
                self_energy[0,:,:,:,0] = self.self_energy[:,:,:].real
                self_energy[0,:,:,:,1] = self.self_energy[:,:,:].imag

            # ZSE fan active
            self_energy_fan_active = ds.createVariable('self_energy_fan_active','d',
                ('number_of_spins', 'number_of_kpoints',
                 'max_number_of_states', 'number_of_frequencies', 'cplex'))

            if self.self_energy_fan_active is not None:

                # FIXME number of spin
                self_energy_fan_active[0,:,:,:,0] = self.self_energy_fan_active[:,:,:].real
                self_energy_fan_active[0,:,:,:,1] = self.self_energy_fan_active[:,:,:].imag

            # ZSE static
            data = ds.createVariable(
                'self_energy_static','d',
                ('number_of_spins', 'number_of_kpoints',
                 'max_number_of_states'))

            if self.self_energy_static is not None:
                # FIXME number of spin
                data[0,:,:] = self.self_energy_static[:,:].real

            # TSE
            self_energy_T = ds.createVariable(
                'self_energy_temperature_dependent','d',
                ('number_of_spins', 'number_of_kpoints',
                 'max_number_of_states', 'number_of_frequencies',
                 'number_of_temperature', 'cplex'))

            if self.self_energy_T is not None:
                # FIXME number of spin
                self_energy_T[0,:,:,:,:,0] = self.self_energy_T[:,:,:,:].real
                self_energy_T[0,:,:,:,:,1] = self.self_energy_T[:,:,:,:].imag

            # TSE static
            data = ds.createVariable(
                'self_energy_static_T','d',
                ('number_of_spins', 'number_of_kpoints',
                 'max_number_of_states', 'number_of_temperature'))

            if self.self_energy_static_T is not None:
                # FIXME number of spin
                data[0,:,:,:] = self.self_energy_static_T[:,:,:].real

            # ZSF
            spectral_function = ds.createVariable(
                'spectral_function','d',
                ('number_of_spins', 'number_of_kpoints',
                 'max_number_of_states', 'number_of_frequencies'))

            if self.spectral_function is not None:
                # FIXME number of spin
                spectral_function[0,:,:,:] = self.spectral_function[:,:,:]

            spectral_function_T = ds.createVariable(
                'spectral_function_temperature_dependent','d',
                ('number_of_spins', 'number_of_kpoints',
                 'max_number_of_states', 'number_of_frequencies',
                 'number_of_temperature'))

            # TSF
            if self.spectral_function_T is not None:
                # FIXME number of spin
                spectral_function_T[0,:,:,:,:] = (
                    self.spectral_function_T[:,:,:,:])
        return


    @master_only
    def write_renormalization(self):
        """Write the computed renormalization in a text file."""
        create_directory(self.ren_dat)

        with open(self.ren_dat, "w") as f:

            if self.zero_point_renormalization is not None:

                f.write("Total zero point renormalization (eV) for "
                        "{} Q points\n".format(self.nqpt))

                for ikpt, kpt in enumerate(self.kpts):
                    f.write('Kpt: {0[0]} {0[1]} {0[2]}\n'.format(kpt))

                    for line in formatted_array_lines(
                        self.zero_point_renormalization[ikpt,:].real*Ha2eV):

                        f.write(line)

            if self.temperature_dependent_renormalization is not None:
                f.write("Temperature dependence at Gamma (eV)\n")

                for iband in range(self.nband):
                  f.write('Band: {}\n'.format(iband))

                  for tt, T in enumerate(self.temperatures):

                    ren = (
                    self.temperature_dependent_renormalization[0,iband,tt]
                    .real * Ha2eV)
                    f.write("{:>8.1f}  {:>12.8f}\n".format(T, ren))

        return

    @master_only
    def write_broadening(self):
        """Write the computed broadening in a text file."""
        create_directory(self.ren_dat)

        with open(self.ren_dat, "w") as f:

            if self.zero_point_broadening is not None:

                f.write("Total zero point broadening (eV) for "
                        "{} Q points\n".format(self.nqpt))

                for ikpt, kpt in enumerate(self.kpts):
                    f.write('Kpt: {0[0]} {0[1]} {0[2]}\n'.format(kpt))
                    for line in formatted_array_lines(
                        self.zero_point_broadening[ikpt,:].real*Ha2eV):

                        f.write(line)

            if self.temperature_dependent_broadening is not None:

                f.write("Temperature dependence at Gamma\n")

                for iband in range(self.nband):
                  f.write('Band: {}\n'.format(iband))

                  for tt, T in enumerate(self.temperatures):

                    brd = (self.temperature_dependent_broadening[0,iband,tt]    
                           .real * Ha2eV)
                    f.write("{:>8.1f}  {:>12.8f}\n".format(T, brd))

        return
