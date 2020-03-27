"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
"HAVE_WANNIER90",
]

#: List of keywords that are automatically added to all the tests of this suite. 
keywords = [
]

#: List of input files
inp_files = [
"t00.in",
"t01.in",
"t02.in",
"t03.in",
"t04.in",
"t11.in",
"t12.in",
"t13.in",
]
