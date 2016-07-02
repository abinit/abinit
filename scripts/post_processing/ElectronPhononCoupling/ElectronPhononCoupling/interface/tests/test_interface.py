
import unittest
import os

import ElectronPhononCoupling

inputsdir = os.path.abspath(os.path.join(__file__, '..', '..', '..', 'data', 'inputs_for_tests'))
datadir = os.path.join(__file__, '..', '..', '..', 'data', 'data_LiF')

class TestCommandLine(unittest.TestCase):

    def setUp(self):
        return

    def tearDown(self):
        os.chdir(inputsdir)
        os.system('rm -r output')

    #def make_input_file(self, ref_input, datadir):
    #    """Read a reference input file, and create a new input with appropriate files."""

    def run_pp_temperature_with_input(self, fname):
        os.chdir(inputsdir)
        os.system('electron-phonon-coupling < ' + fname)

    #def test_test11(self): self.run_pp_temperature_with_input('t11.in')
    #def test_test12(self): self.run_pp_temperature_with_input('t12.in')
    #def test_test13(self): self.run_pp_temperature_with_input('t13.in')
    #def test_test14(self): self.run_pp_temperature_with_input('t14.in')
    #def test_test21(self): self.run_pp_temperature_with_input('t21.in')
    #def test_test22(self): self.run_pp_temperature_with_input('t22.in')
    #def test_test23(self): self.run_pp_temperature_with_input('t23.in')
    #def test_test24(self): self.run_pp_temperature_with_input('t24.in')
    #def test_test31(self): self.run_pp_temperature_with_input('t31.in')
    #def test_test32(self): self.run_pp_temperature_with_input('t32.in')
    #def test_test33(self): self.run_pp_temperature_with_input('t33.in')
    #def test_test34(self): self.run_pp_temperature_with_input('t34.in')
    #def test_test41(self): self.run_pp_temperature_with_input('t41.in')
