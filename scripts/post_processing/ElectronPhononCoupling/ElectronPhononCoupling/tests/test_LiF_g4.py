from __future__ import print_function
from os.path import join as pjoin
from copy import copy

from . import SETest

from ..data.LiF_g4 import nqpt, wtq, fnames, refdir

class Test_LiF_g4(SETest):

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

        nqpt = nqpt,
        wtq = wtq,
        smearing_eV = 0.01,
        temp_range = [0, 300, 300],
        omega_range = [-0.1, 0.1, 0.001],
        rootname = 'epc.out',
        **fnames)

    @property
    def refdir(self):
        return refdir

    # ZPR
    def test_zpr_dyn(self):
        """Dynamical zero-point renormalization"""
        self.run_compare_nc(
            function = self.get_zpr_dyn,
            key = 'zero_point_renormalization',
            )

    #def generate_zpr_dyn(self):
    #    """Generate epc data for this test."""
    #    return self.generate_test_ref(self.get_zpr_dyn)

    def test_tdr_dyn(self):
        """Dynamical temperature dependent renormalization"""
        self.run_compare_nc(
            function = self.get_tdr_dyn,
            key = 'temperature_dependent_renormalization',
            )

    # ZPR
    def test_zpr_stat_mode(self):
        """Dynamical zero-point renormalization"""
        self.run_compare_nc(
            function = self.get_zpr_stat_mode,
            key = 'zero_point_renormalization_by_modes',
            )

    def test_zpb_dyn(self):
        """Dynamical ZP Brd"""
        self.run_compare_nc(
            function = self.get_zpb_dyn,
            key = 'zero_point_broadening',
            )

    def test_tdb_dyn(self):
        """Dynamical TD Brd"""
        self.run_compare_nc(
            function = self.get_tdb_dyn,
            key = 'temperature_dependent_broadening',
            )

    def test_zpb_stat(self):
        """Static ZP Brd"""
        self.run_compare_nc(
            function = self.get_zpb_stat,
            key = 'zero_point_broadening',
            )

    def test_tdb_stat(self):
        """Dynamical TD Brd"""
        self.run_compare_nc(
            function = self.get_tdb_stat,
            key = 'temperature_dependent_broadening',
            )

    #def generate_tdr_dyn(self):
    #    return self.generate_test_ref(self.get_tdr_dyn)

    # All
    def generate(self):
        """Generate epc data for all tests."""

        print('Generating reference data for tests in directory: {}'.format(
              self.refdir))

        for function in (
            self.get_zpr_dyn,
            self.get_tdr_dyn,
            self.get_zpr_stat_mode,
            self.get_zpb_dyn,
            self.get_tdb_dyn,
            self.get_zpb_stat,
            self.get_tdb_stat,
            ):
            self.generate_ref(function)


