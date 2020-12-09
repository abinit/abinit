"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
"!HAVE_MPI",
]

#: List of keywords that are automatically added to all the tests of this suite. 
keywords = [
]

subsuites = [
"sv2",
"sv3",
"sv4",
"sv5",
"sv6",
"sv7",
]


#: List of input files
inp_files = [
"tsv2_81.abi",
"tsv2_82.abi",
"tsv3_03.abi",
"tsv3_04.abi",
"tsv3_05.abi",
"tsv4_55.abi",
"tsv4_78.abi",
"tsv4_80.abi",
"tsv4_90.abi",
"tsv5_112.abi",
"tsv5_113.abi",
"tsv6_121.abi",
"tsv6_122.abi",
"tsv6_123.abi",
"tsv6_124.abi",
"tsv6_125.abi",
"tsv6_126.abi",
"tsv7_70.abi",
]
