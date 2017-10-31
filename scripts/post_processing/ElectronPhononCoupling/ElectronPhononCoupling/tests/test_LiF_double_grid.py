from __future__ import print_function
from os.path import join as pjoin
from copy import copy

from . import SETest

from ..data import LiF_g4 as g4
from ..data import LiF_g2_2 as g2

class Test_LiF_double_grid(SETest):

    common = dict(
        temperature = False,
        renormalization = False,
        broadening = False,
        self_energy = False,
        spectral_function = False,
        dynamical = True,
        split_active = True,
        double_grid = True,
        write = True,
        verbose = False,

        nqpt = g2.nqpt,
        wtq = g2.wtq,

        nqpt_fine = g4.nqpt,
        wtq_fine = g4.wtq,

        smearing_eV = 0.01,
        temp_range = [0, 300, 300],
        omega_range = [-0.1, 0.1, 0.001],

        eigq_fine_fnames = g4.fnames['eigq_fnames'],
        ddb_fine_fnames = g4.fnames['ddb_fnames'],
        gkk_fine_fnames = g4.fnames['gkk_fnames'],
        **g2.fnames)

    @property
    def refdir(self):
        return g2.refdir

    def test_tdr_dyn_double_grid(self):
        """Dynamical temperature dependent renormalization"""
        self.run_compare_nc(
            function = self.get_tdr_dyn_double_grid,
            key = 'temperature_dependent_renormalization',
            )

    def test_zp_se_double_grid(self):
        """Dynamical temperature dependent renormalization"""
        self.run_compare_nc(
            function = self.get_zp_se_double_grid,
            key = 'self_energy',
            )

    def test_td_se_double_grid(self):
        """Dynamical temperature dependent renormalization"""
        self.run_compare_nc(
            function = self.get_td_se_double_grid,
            key = 'self_energy_temperature_dependent',
            )

    # All
    def generate(self):
        """Generate epc data for all tests."""

        print('Generating reference data for tests in directory: {}'.format(
              self.refdir))

        for function in (
            self.get_tdr_dyn_double_grid,
            self.get_zp_se_double_grid,
            self.get_td_se_double_grid,
            ):
            self.generate_ref(function)


