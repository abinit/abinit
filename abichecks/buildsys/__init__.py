"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
]

#: List of keywords that are automatically added to all the tests of this suite. 
keywords = [
"buildsys",
]

#: List of input files
#inp_files = [
#]

#: List of python scripts. 
pyscripts = [
"check-binaries-conf.py",
"check-build-config.py",
"check-cpp-options.py",
"check-forbidden-flags.py",
]
