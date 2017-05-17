
from . import test_LiF_g2
from . import test_LiF_g4
from . import test_LiF_double_grid

def generate_all_data_for_tests():
    """Generate the epc data for all tests."""

    print('Generating all data for the tests suite.')
    test_LiF_g2.Test_LiF_g2('generate').generate()
    test_LiF_g4.Test_LiF_g4('generate').generate()
    test_LiF_double_grid.Test_LiF_double_grid('generate').generate()


