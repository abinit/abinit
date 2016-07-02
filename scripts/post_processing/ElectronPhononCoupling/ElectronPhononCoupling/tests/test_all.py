from os.path import join as pjoin
from copy import copy

from . import EPCTest
from ..data import LiF_nqpt, LiF_wtq, LiF_fnames, LiF_outputs

from .. import compute_epc


class AllTests(EPCTest):

    common = dict(
        write=True,
        smearing_eV=0.01,
        temp_range=[0,600,300],
        nqpt=LiF_nqpt,
        wtq=LiF_wtq,
        **LiF_fnames)

    def test_t11(self):
        """Static ZP Ren"""

        root = pjoin(self.tmpdir, 't11')
        ref = LiF_outputs['t11']
        out = root + '_EP.nc'

        compute_epc(
            calc_type=1,
            temperature=False,
            lifetime=False,
            output=root,
            **self.common)

        self.AssertClose(out, ref, 'zero_point_renormalization')

    def test_t12(self):
        """Static Tdep Ren"""

        root = pjoin(self.tmpdir, 't12')
        ref = LiF_outputs['t12']
        out = root + '_EP.nc'

        compute_epc(
            calc_type=1,
            temperature=True,
            lifetime=False,
            output=root,
            **self.common)

        self.AssertClose(out, ref, 'temperature_dependent_renormalization')

    def test_t13(self):
        """Static ZP Brd"""

        root = pjoin(self.tmpdir, 't13')
        ref = LiF_outputs['t13']
        out = root + '_EP.nc'

        compute_epc(
            calc_type=1,
            temperature=False,
            lifetime=True,
            output=root,
            **self.common)

        self.AssertClose(out, ref, 'zero_point_broadening')

    def test_t14(self):
        """Static Tdep Brd"""

        root = pjoin(self.tmpdir, 't14')
        ref = LiF_outputs['t14']
        out = root + '_EP.nc'

        compute_epc(
            calc_type=1,
            temperature=True,
            lifetime=True,
            output=root,
            **self.common)

        self.AssertClose(out, ref, 'temperature_dependent_broadening')

    def test_t21(self):
        """Static ZP Ren"""

        root = pjoin(self.tmpdir, 't21')
        ref = LiF_outputs['t21']
        out = root + '_EP.nc'

        compute_epc(
            calc_type=2,
            temperature=False,
            lifetime=False,
            output=root,
            **self.common)

        self.AssertClose(out, ref, 'zero_point_renormalization')

    def test_t22(self):
        """Static Tdep Ren"""

        root = pjoin(self.tmpdir, 't22')
        ref = LiF_outputs['t22']
        out = root + '_EP.nc'

        compute_epc(
            calc_type=2,
            temperature=True,
            lifetime=False,
            output=root,
            **self.common)

        self.AssertClose(out, ref, 'temperature_dependent_renormalization')

    def test_t23(self):
        """Static ZP Brd"""

        root = pjoin(self.tmpdir, 't23')
        ref = LiF_outputs['t23']
        out = root + '_EP.nc'

        compute_epc(
            calc_type=2,
            temperature=False,
            lifetime=True,
            output=root,
            **self.common)

        self.AssertClose(out, ref, 'zero_point_broadening')

    def test_t24(self):
        """Static Tdep Brd"""

        root = pjoin(self.tmpdir, 't24')
        ref = LiF_outputs['t24']
        out = root + '_EP.nc'

        compute_epc(
            calc_type=2,
            temperature=True,
            lifetime=True,
            output=root,
            **self.common)

        self.AssertClose(out, ref, 'temperature_dependent_broadening')

    def test_t31(self):
        """Static ZP Ren"""

        root = pjoin(self.tmpdir, 't31')
        ref = LiF_outputs['t31']
        out = root + '_EP.nc'

        compute_epc(
            calc_type=3,
            temperature=False,
            lifetime=False,
            output=root,
            **self.common)

        self.AssertClose(out, ref, 'zero_point_renormalization')

    def test_t32(self):
        """Static Tdep Ren"""

        root = pjoin(self.tmpdir, 't32')
        ref = LiF_outputs['t32']
        out = root + '_EP.nc'

        compute_epc(
            calc_type=3,
            temperature=True,
            lifetime=False,
            output=root,
            **self.common)

        self.AssertClose(out, ref, 'temperature_dependent_renormalization')

    def test_t33(self):
        """Static ZP Brd"""

        root = pjoin(self.tmpdir, 't33')
        ref = LiF_outputs['t33']
        out = root + '_EP.nc'

        compute_epc(
            calc_type=3,
            temperature=False,
            lifetime=True,
            output=root,
            **self.common)

        self.AssertClose(out, ref, 'zero_point_broadening')

    def test_t34(self):
        """Static Tdep Brd"""

        root = pjoin(self.tmpdir, 't34')
        ref = LiF_outputs['t34']
        out = root + '_EP.nc'

        compute_epc(
            calc_type=3,
            temperature=True,
            lifetime=True,
            output=root,
            **self.common)

        self.AssertClose(out, ref, 'temperature_dependent_broadening')

