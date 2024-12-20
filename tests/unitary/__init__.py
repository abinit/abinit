"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
]

#: List of keywords that are automatically added to all the tests of this suite. 
keywords = [
]

subsuites = [
"fftsg",
"fftfftw3",
"fftmkl",
"fftgw",
"fourdp",
"fourwf",
"transposer",
]

#: List of input files
inp_files = [
"tfftsg_01.abi",
"tfftsg_02.abi",
"tfftsg_03.abi",
"tfftsg_04.abi",
"tfftsg_05.abi",
"tfftsg_06.abi",
#
"tfftfftw3_01.abi",
"tfftfftw3_02.abi",
"tfftfftw3_03.abi",
"tfftfftw3_04.abi",
"tfftfftw3_05.abi",
"tfftfftw3_06.abi",
#
"tfftmkl_01.abi",
"tfftmkl_02.abi",
"tfftmkl_03.abi",
"tfftmkl_04.abi",
#
"tfftgw_01.abi",
"tfftgw_02.abi",
"tfftgw_03.abi",
#
"tfourdp_01.abi",
"tfourdp_02.abi",
#
"tfourwf_01.abi",
"tfourwf_02.abi",
"tfourwf_03.abi",
"tfourwf_04.abi",
"tfourwf_05.abi",
"tfourwf_06.abi",
"tfourwf_07.abi",
"tfourwf_08.abi",
"tfourwf_09.abi",
"tfourwf_10.abi",
"tfourwf_11.abi",
"tfourwf_12.abi",
"tfourwf_13.abi",
#
"ttransposer_01.abi",
]
