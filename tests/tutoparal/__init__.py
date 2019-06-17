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
#"gswvl_01.in",
#"gswvl_02.in",
#"tdfpt_01.in",
#"tdfpt_02.in",
"tdfpt_03.in",
"tdfpt_04.in",
#"tdfpt_03PAW.in,"
#"tdfpt_04PAW.in",
#"tdmft_1.in",
#"tdmft_2.in",
"tgspw_01.in",
"tgspw_02.in",    # OK
"tgspw_03.in",    # OK
#"tgspw_04.in",   # Unstable because nstep=5 and bandpp: 2d iteration oscillates and fldiff does not handle it!
#"tgspw_05.in",
"tmbt_1.in",   # OK     
"tmbt_2.in",   # OK
"tmbt_3.in",   # OK
"tmbt_4.in",   # OK
#"tmbt_5.in",  
#"tmoldyn_01.in",    # This is not stable. Does it use random number generators?
#"tmoldyn_02.in",
#"tmoldyn_03.in",
#"tmoldyn_04.in",
#"tmoldyn_05.in",
#"tmoldyn_06.in",
#"tmoldyn_07.in",
"tstring_01.in",   # MPI Error in MPI_ALLREDUCE
#"tstring_02.in",
#"tstring_03.in",
#"tstring_04.in",
"tucrpa_1.in",
"tucrpa_2.in",
#"tucrpa_3.in",
#"tucrpa_4.in",
#"tucrpa_5.in",
]
