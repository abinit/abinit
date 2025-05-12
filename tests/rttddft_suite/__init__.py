"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
]

#: List of keywords that are automatically added to all the tests of this suite.
keywords = [
"RTTDDFT",
]

#: List of input files
inp_files = [
"t01.abi", # Tests ER and EMR propagator in PAW (Ni FCC)
"t02.abi", # Tests ER and EMR propagator in NC (Al FCC)
"t03.abi", # Tests restart capabilities in stationary case (Al FCC NC)
"t04.abi", # Tests restart capabilities in non stationary case (Ni FCC PAW)
"t05.abi", # Tests Impulsional TD electric field without induced potential vector (C Diamond 1 atom PAW)
"t06.abi", # Tests Impulsional TD electric field with induced potential vector (C Diamond 1 atom PAW)
]
