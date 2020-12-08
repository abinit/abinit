"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
"HAVE_MPI",
]

#: List of keywords that are automatically added to all the tests of this suite. 
keywords = [
]

#: This suite contains tests executed with different numbers of MPI processes.
is_multi_parallel = True

subsuites = [
"dfpt",
"gspw",
#"gswvl",
"mbt",
#"dmft",
#"moldyn",
"string",
"ucrpa",
]

#: List of input files
inp_files = [
#"gswvl_01.abi",
#"gswvl_02.abi",
#"tdfpt_01.abi",
#"tdfpt_02.abi",
"tdfpt_03.abi",
"tdfpt_04.abi",
#"tdfpt_03PAW.abi,"
#"tdfpt_04PAW.abi",
#"tdmft_1.abi",
#"tdmft_2.abi",
"tgspw_01.abi",
"tgspw_02.abi",    # OK
"tgspw_03.abi",    # OK
#"tgspw_04.abi",   # Unstable because nstep=5 and bandpp: 2d iteration oscillates and fldiff does not handle it!
#"tgspw_05.abi",
"tmbt_1.abi",   # OK     
"tmbt_2.abi",   # OK
"tmbt_3.abi",   # OK
"tmbt_4.abi",   # OK
#"tmbt_5.abi",  
#"tmoldyn_01.abi",    # This is not stable. Does it use random number generators?
#"tmoldyn_02.abi",
#"tmoldyn_03.abi",
#"tmoldyn_04.abi",
#"tmoldyn_05.abi",
#"tmoldyn_06.abi",
#"tmoldyn_07.abi",
"tstring_01.abi",   # MPI Error in MPI_ALLREDUCE
#"tstring_02.abi",
#"tstring_03.abi",
#"tstring_04.abi",
"tucrpa_1.abi",
"tucrpa_2.abi",
#"tucrpa_3.abi",
#"tucrpa_4.abi",
#"tucrpa_5.abi",
]
