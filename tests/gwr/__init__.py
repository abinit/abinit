"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
"HAVE_LINALG_SCALAPACK",
]

#: List of keywords that are automatically added to all the tests of this suite.
keywords = [
"GWR",
]

#: List of input files
inp_files = [
"-t01.abi", # disabled for now
"-t02.abi", # disabled for now
"-t03.abi", # disabled for now
"-t04.abi", # disabled for now
"-t05.abi", # disabled for now
"-t06.abi", # disabled for now
"-t07.abi", # disabled for now
]
