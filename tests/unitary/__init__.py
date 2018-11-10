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
"tfftsg_01.in",
"tfftsg_02.in",
"tfftsg_03.in",
"tfftsg_04.in",
"tfftsg_05.in",
"tfftsg_06.in",
#
"tfftfftw3_01.in",
"tfftfftw3_02.in",
"tfftfftw3_03.in",
"tfftfftw3_04.in",
"tfftfftw3_05.in",
"tfftfftw3_06.in",
#
"tfftmkl_01.in",
"tfftmkl_02.in",
"tfftmkl_03.in",
"tfftmkl_04.in",
#
"tfftgw_01.in",
"tfftgw_02.in",
"tfftgw_03.in",
#
"tfourdp_01.in",
"tfourdp_02.in",
#
"tfourwf_01.in",
"tfourwf_02.in",
"tfourwf_03.in",
"tfourwf_04.in",
"tfourwf_05.in",
"tfourwf_06.in",
"tfourwf_07.in",
"tfourwf_08.in",
"tfourwf_09.in",
"tfourwf_10.in",
"tfourwf_11.in",
"tfourwf_12.in",
"tfourwf_13.in",
#
"ttransposer_01.in",
]
