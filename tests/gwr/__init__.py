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
"t01.abi",
"t02.abi",
"t03.abi",
"t04.abi",
"t05.abi",
"t06.abi",
"t07.abi",
"t08.abi",
"t09.abi",
"t10.abi" ,  # test Gamma-only HF with time-reversal symmetry (real wavefunctions)
"t11.abi" ,  # test Gamma-only HF without time-reversal symmetry (complex wavefunctions)
]
