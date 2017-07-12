from os.path import join as pjoin
from copy import copy

from . import EPCTest, SETest
from ..data import LiF_g2 as test


# FIXME

class Test_LiF_g2(SETest):

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

        nqpt=test.nqpt,
        wtq=test.wtq,
        smearing_eV=0.01,
        temp_range=[0,300,300],
        omega_range=[-.1,.1,.001],

        rootname = 'epc.out',
        **test.fnames)

    @property
    def refdir(self):
        return test.refdir

    def test_zpr_dyn(self):
        """Dynamical ZPR"""
        self.run_compare_nc(
            function = self.get_zpr_dyn,
            key = 'zero_point_renormalization',
            )

    def test_tdr_dyn(self):
        """Dynamical Tdep Ren"""
        self.run_compare_nc(
            function = self.get_tdr_dyn,
            key = 'temperature_dependent_renormalization',
            )

    def test_zp_se(self):
        """Zero Point Self-Energy"""
        self.run_compare_nc(
            function = self.get_zp_se,
            key = 'self_energy',
            )

    def test_zp_sf(self):
        """Zero Point Spectral Function"""
        self.run_compare_nc(
            function = self.get_zp_sf,
            key = 'spectral_function',
            )

    def test_td_se(self):
        """Temperature Dependent Self-Energy"""
        self.run_compare_nc(
            function = self.get_td_se,
            key = 'self_energy_temperature_dependent',
            )

    def test_td_sf(self):
        """Temperature Dependent Spectral Function"""
        self.run_compare_nc(
            function = self.get_td_sf,
            key = 'spectral_function_temperature_dependent',
            )

    def test_zpr_stat(self):
        """Static ZP Ren"""
        self.run_compare_nc(
            function = self.get_zpr_stat,
            key = 'zero_point_renormalization',
            )

    def test_tdr_stat(self):
        """Static Tdep Ren"""
        self.run_compare_nc(
            function = self.get_tdr_stat,
            key = 'temperature_dependent_renormalization',
            )

    def test_zpr_stat_nosplit(self):
        """Static Zero Point Renormalization"""
        self.run_compare_nc(
            function = self.get_zpr_stat_nosplit,
            key = 'zero_point_renormalization',
            )

    def test_tdr_static_nosplit(self):
        """Static Temperature Dependent Renormalization"""
        self.run_compare_nc(
            function = self.get_tdr_stat_nosplit,
            key = 'temperature_dependent_renormalization',
            )

    def test_zpb_stat_nosplit(self):
        """Static Zero Point Broadening"""
        self.run_compare_nc(
            function = self.get_zpb_stat_nosplit,
            key = 'zero_point_broadening',
            )

    def test_tdb_stat_nosplit(self):
        """Static Temperature Dependent Broadening"""
        self.run_compare_nc(
            function = self.get_tdb_stat_nosplit,
            key = 'temperature_dependent_broadening',
            )

    # All
    def generate(self):
        """Generate epc data for all tests."""

        print('Generating reference data for tests in directory: {}'.format(
              self.refdir))

        for function in (
            self.get_zpr_dyn,
            self.get_tdr_dyn,
            self.get_zp_se,
            self.get_zp_sf,
            self.get_td_se,
            self.get_td_sf,
            self.get_zpr_stat,
            self.get_tdr_stat,
            self.get_zpr_stat_nosplit,
            self.get_tdr_stat_nosplit,
            self.get_zpb_stat_nosplit,
            self.get_tdb_stat_nosplit,
            ):
            self.generate_ref(function)


