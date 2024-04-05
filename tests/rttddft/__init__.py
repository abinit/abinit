"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
]

#: List of keywords that are automatically added to all the tests of this suite.
keywords = [
]

#: List of input files
inp_files = [
"t01.abi", # Tests ER and EMR propagator in PAW (Ni FCC nonmagnetic)
"t02.abi", # Tests ER and EMR propagator in NC (Al FCC)
"t03.abi", # Tests restart capabilities (ER propagator Al FCC NC)
]
