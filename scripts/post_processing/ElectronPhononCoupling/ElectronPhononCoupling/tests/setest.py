
import os
from os.path import join as pjoin
from copy import copy

from ..core.mpi import master_only
from ..interface import compute

from .epctest import EPCTest

__all__ = ['SETest']

class SETest(EPCTest):
    """
    Base class for tests involving the electron self-energy.
    """

    common = dict(
        temperature = False,
        renormalization = False,
        broadening = False,
        self_energy = False,
        spectral_function = False,
        dynamical = True,
        split_active = True,
        double_grid = False,
        write = True,
        verbose = False,

        nqpt = 1,
        wtq = [1.],
        smearing_eV = 0.01,
        temp_range = [0, 300, 300],
        omega_range = [-0.1, 0.1, 0.001],
        rootname = 'epc.out',
        )

    def get_kwargs(self, dirname, basename, **kwargs):
        """Construct the input dictionary for a test"""
        new_kwargs = copy(self.common)
        new_kwargs.update(rootname=pjoin(dirname, basename))
        new_kwargs.update(**kwargs)
        return new_kwargs

    def run_compare_nc(self, function, key, refdir=None, nc_ref=None):
        """
        Run 'compute' by generating the arguments with 'function'
        then compare the array 'key' in nc_output.

        key:
            Key to be compared in both nc_output.
        function:
            Function generating the keyword arguments for 'compute'.
            Takes a directory as single argument.
        refdir:
            directory for reference files
        nc_ref:
            name of the netcdf file for comparison.
            Alternate argument of refdir.
        """

        out = compute(**function(self.tmpdir))

        if nc_ref is not None:
            pass
        else:
            if refdir is None:
                refdir = self.refdir
            nc_ref = out.nc_output.replace(self.tmpdir, refdir)

        self.check_reference_exists(nc_ref)
        self.AssertClose(out.nc_output, nc_ref, key)

    def generate_ref(self, function):
        """
        function:
            Function generating the keyword arguments for 'compute'.
            Takes a directory as single argument.
        """
        return compute(**function(self.refdir))

    def get_zpr_dyn(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='zpr_dyn',
            renormalization=True,
            )

    def get_tdr_dyn(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='tdr_dyn',
            renormalization=True,
            temperature=True,
            )

    def get_zpb_dyn(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='zpb_dyn',
            temperature=False,
            broadening=True,
            dynamical=True,
            )

    def get_zp_se(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='zp_se',
            temperature=False,
            self_energy=True,
            )

    def get_td_se(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='td_se',
            temperature=True,
            self_energy=True,
            )

    def get_zp_sf(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='zp_sf',
            temperature=False,
            self_energy=True,
            spectral_function=True,
            )

    def get_td_sf(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='td_sf',
            temperature=True,
            self_energy=True,
            spectral_function=True,
            )


    def get_zpr_stat(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='zpr_stat',
            temperature=False,
            renormalization=True,
            dynamical=False,
            )

    def get_zpr_stat_mode(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='zpr_stat_mode',
            temperature=False,
            renormalization=True,
            dynamical=False,
            mode=True,
            )

    def get_tdr_stat(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='tdr_stat',
            temperature=True,
            renormalization=True,
            dynamical=False,
            )

    def get_zpb_stat(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='zpb_stat',
            temperature=False,
            broadening=True,
            dynamical=False,
            )

    def get_tdb_stat(self, dirname):
        return self.get_kwargs(
            basename='tdb_stat',
            temperature=True,
            broadening=True,
            dynamical=False,
            )

    def get_zpr_stat_nosplit(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='zpr_stat_nosplit',
            temperature=False,
            renormalization=True,
            dynamical=False,
            split_active=False,
            )

    def get_tdr_stat_nosplit(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='tdr_stat_nosplit',
            temperature=True,
            renormalization=True,
            dynamical=False,
            split_active=False,
            )

    def get_zpb_dyn(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='zpb_dyn',
            temperature=False,
            broadening=True,
            dynamical=True,
            split_active=True,
            )

    def get_tdb_dyn(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='tdb_dyn',
            temperature=True,
            broadening=True,
            dynamical=True,
            split_active=True,
            )

    def get_zpb_stat(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='zpb_stat',
            temperature=False,
            broadening=True,
            dynamical=False,
            split_active=True,
            )

    def get_tdb_stat(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='tdb_stat',
            temperature=True,
            broadening=True,
            dynamical=False,
            split_active=True,
            )

    def get_zpb_stat_nosplit(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='zpb_stat_nosplit',
            temperature=False,
            broadening=True,
            dynamical=False,
            split_active=False,
            )

    def get_tdb_stat_nosplit(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='tdb_stat_nosplit',
            temperature=True,
            broadening=True,
            dynamical=False,
            split_active=False,
            )

    def get_tdr_dyn_double_grid(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='tdb_dyn_double_grid',
            temperature=True,
            renormalization=True,
            dynamical=True,
            split_active=True,
            double_grid=True,
            )

    def get_zp_se_double_grid(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='zp_se_double_grid',
            self_energy=True,
            temperature=False,
            double_grid=True,
            )

    def get_td_se_double_grid(self, dirname):
        return self.get_kwargs(
            dirname,
            basename='td_se_double_grid',
            self_energy=True,
            temperature=True,
            double_grid=True,
            )

