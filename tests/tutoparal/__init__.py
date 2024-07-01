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
"dmft",
"paral_bandpw",
"gswvl",
"images",
"mbt",
"moldyn",
"psic",
"ucalc_crpa",
]

#: List of input files
inp_files = [
"tdfpt_01.abi",
"tdfpt_02.abi",
"tdfpt_03.abi",
"tdfpt_04.abi",
#"tdfpt_03PAW.abi,"
#"tdfpt_04PAW.abi",
"tdmft_1.abi",
"tdmft_2.abi",
"tparal_bandpw_01.abi",
"tparal_bandpw_02.abi",    # OK
"tparal_bandpw_03.abi",    # OK
"tparal_bandpw_04.abi",   # Unstable because nstep=5 and bandpp: 2d iteration oscillates and fldiff does not handle it!
#"tparal_bandpw_05.abi",
"tgswvl_1.abi",
"tgswvl_2.abi",
"timages_01.abi",
"timages_02.abi",
"timages_03.abi",
"timages_04.abi",
"tmbt_1.abi",   # OK     
"tmbt_2.abi",   # OK
"tmbt_3.abi",   # OK
#"tmbt_4.abi",   # OK on some machines, but not on S64 .
#"tmbt_5.abi",  
 "tmoldyn_01.abi",    # Not really tested on all slaves, runs are a bit long. 
#"tmoldyn_02.abi",
#"tmoldyn_03.abi",
#"tmoldyn_04.abi",
#"tmoldyn_05.abi",    # Run is too long : more than 20 minutes on 64 procs
#"tmoldyn_06.abi",
#"tmoldyn_07.abi",    # Run is too long : more than 20 minutes on 2 procs.
"tpsic_01.abi", 
"tpsic_02.abi",   
"tpsic_03.abi", 
"tucalc_crpa_1.abi",
"tucalc_crpa_2.abi",
#"tucalc_crpa_3.abi",
#"tucalc_crpa_4.abi",
#"tucalc_crpa_5.abi",
]
