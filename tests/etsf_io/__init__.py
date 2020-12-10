"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
"HAVE_NETCDF",
]

#: List of keywords that are automatically added to all the tests of this suite.
keywords = [
"netcdf"
]

#: List of input files
inp_files = [
"t00.abi",
"t01.abi",
"t02.abi",
"t03.abi",
"t04.abi",
"t09.abi",
"-t10.abi", # Disabled
"t21.abi",
"t22.abi",
"t30.abi",
]
