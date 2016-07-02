"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
]

#: List of keywords that are automatically added to all the tests of this suite. 
keywords = [
"abirules",
]

#: List of input files
#inp_files = [
#]

#: List of python scripts. 
pyscripts = [
###"check-input-vars.py",
"check_ascii.py",
"check_forbidden.py",
"check_forbidden_in_doc_yml.py",
"check_inlined_macros.py",
"docchk.py",
"warningschk.py",
]
