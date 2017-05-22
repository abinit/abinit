import os
from os.path import join as pjoin
from copy import copy
import unittest
import tempfile, shutil

import numpy as np 
from numpy.testing import assert_allclose
import netCDF4 as nc

from ..core.mpi import master_only

__all__ = ['EPCTest']

class EPCTest(unittest.TestCase):
    """Base class for tests."""

    #def runTest(self):
    #    """Dummy method."""
    #    pass

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.local_testdir = 'tmp.Tests'
        self.local_tmpdir = os.path.join(
            self.local_testdir, os.path.split(self.tmpdir)[-1])

    @property
    def refdir(self):
        raise NotImplementedError()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def recover_tmpdir(self):
        """Debug tool to copy the execution directory locally."""
        if not os.path.exists(self.local_testdir):
            os.mkdir(self.local_testdir)
        shutil.copytree(self.tmpdir, self.local_tmpdir)

    #@master_only  # This decorator screws up the output of nosetests...
    def AssertClose(self, f1, f2, key, **kwargs):
        """Assert that an array in two nc files is close."""

        with nc.Dataset(f1, 'r') as ds1:
            a1 = ds1.variables[key][...]

        with nc.Dataset(f2, 'r') as ds2:
            a2 = ds2.variables[key][...]

        kwargs.setdefault('rtol', 1e-5)
        return assert_allclose(
            a1, a2, err_msg=('\n' +
                '-----------------------\n' +
                'Comparison test failed.\n' +
                'f1: {}\n'.format(f1) +
                'f2: {}\n'.format(f2) +
                'Variable compared: {}\n'.format(key) +
                '-----------------------\n'
                ),
            **kwargs)

    def check_reference_exists(self, fname):
        if not os.path.exists(fname):
            raise Exception('Reference file does not exist: {}'.format(fname))

