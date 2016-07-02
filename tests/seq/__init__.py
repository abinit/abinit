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
"tsv2_81.in",
"tsv2_82.in",
"tsv3_03.in",
"tsv3_04.in",
"tsv3_05.in",
"tsv4_55.in",
"tsv4_78.in",
"tsv4_80.in",
"tsv4_90.in",
"tsv5_112.in",
"tsv5_113.in",
"tsv6_121.in",
"tsv6_122.in",
"tsv6_123.in",
"tsv6_124.in",
"tsv6_125.in",
"tsv6_126.in",
"tsv7_70.in",
]
