"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
]

#: List of keywords that are automatically added to all the tests of this suite. 
keywords = [
]

subsuites = [
"elast",
"depes",
"eph",
"eph4mob",
"ffield",
"lw",
"nlo",
"optic",   
"rf1",
"rf2",
]

#: List of input files
inp_files = [
"tdepes_1.in",
"tdepes_2.in",
"tdepes_3.in",
"tdepes_4.in",
"telast_1.in",
"telast_2.in",
"telast_3.in",
"telast_4.in",
"telast_5.in",
"telast_6.in",
"teph_1.in",
"teph_2.in",
"teph_3.in",
"teph_4.in",
"teph_5.in",
"teph_6.in",
"teph4mob_1.in",
"teph4mob_2.in",
"teph4mob_3.in",
"teph4mob_4.in",
"teph4mob_5.in",
"teph4mob_6.in",
"teph4mob_7.in",
"tffield_1.in",
"tffield_2.in",
"tffield_3.in",
"tffield_4.in",
"tffield_5.in",
"tffield_6.in",
#"-tffield_7.in",  # Disabled
"tlw_1.in",
"tlw_2.in",
"tlw_3.in",
"tlw_4.in",
"tlw_5.in",
"tlw_6.in",
"tlw_7.in",
"tnlo_1.in",
"tnlo_2.in",
"tnlo_3.in",
"tnlo_4.in",
"tnlo_5.in",
"tnlo_6.in",      
"tnlo_7.in",
"tnlo_8.in",
"tnlo_9.in",
"tnlo_10.in",
"tnlo_11.in",
"toptic_1.in",   
"toptic_2.in",
"toptic_3.in",
"toptic_4.in",
"toptic_5.in",
"trf1_1.in",
"trf1_2.in",
"trf1_3.in",
"trf1_4.in",
"trf1_5.in",
"trf1_6.in",
"trf2_1.in",
"trf2_2.in",
"trf2_3.in",
"trf2_4.in",
"trf2_5.in",   # FIXME recheck band2eps
"trf2_6.in",
"trf2_7.in",
]
